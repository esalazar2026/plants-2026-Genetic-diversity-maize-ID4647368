# =============================================================================
# Core collection ESTRATIFICADA con Core Hunter 3 (paquete R 'corehunter')
# Maiz DArTseq - 1,157 accesiones, 20,956 SNPs
#
# Archivos de entrada (en la misma carpeta que este script):
#   genotypes_biparental_20956SNP_1157acc.csv.gz  -> 0/1/2 (1 = heterocigoto), faltantes = NA
#   strata_core_1157acc.csv                       -> estrato y representantes de raza
#
# Diseno:
#   - 2 estratos: LANDRACE (978) y BRED MATERIAL (179). Los 47 UNKNOWN ya vienen
#     asignados segun el tipo inferido por k-NN (40 landrace, 7 bred material).
#   - Cupos proporcionales: 195 landrace + 36 bred material = 231 (20%).
#   - 22 representantes de raza (1 por raza) forzados con always.selected dentro
#     del estrato LANDRACE (cuentan dentro de los 195).
#   - Todas las accesiones son elegibles (sin exclusion por datos faltantes).
#   - Se corren varias configuraciones de objetivos; cada una produce un CSV.
# =============================================================================

# install.packages("corehunter")   # requiere Java 8 o superior
library(corehunter)

# ------------------------------ PARAMETROS ----------------------------------
SEED       <- 20260928           # semilla fija (reproducibilidad)
SIZE       <- c("LANDRACE" = 195, "BRED MATERIAL" = 36)
IMPR_STEPS <- 2000               # termina tras N pasos sin mejora (criterio por pasos, no por tiempo)

# Configuraciones a probar: nombre = lista de objetivos
CONFIGS <- list(
  CV07_EN03      = list(objective("CV", weight = 0.7), objective("EN", "MR", weight = 0.3)),
  CV03_AN07      = list(objective("CV", weight = 0.3), objective("AN", "MR", weight = 0.7)),
  CV02_AN04_EN04 = list(objective("CV", weight = 0.2), objective("AN", "MR", weight = 0.4),
                        objective("EN", "MR", weight = 0.4)),
  CV05_AN05      = list(objective("CV", weight = 0.5), objective("AN", "MR", weight = 0.5))
)

# ------------------------------ DATOS ---------------------------------------
geno_all <- read.csv(gzfile("genotypes_biparental_20956SNP_1157acc.csv.gz"),
                     row.names = 1, check.names = FALSE, na.strings = "NA")
strata   <- read.csv("strata_core_1157acc.csv", stringsAsFactors = FALSE)
strata$ID <- trimws(as.character(strata$ID))
# Conversion explicita a logico (acepta TRUE/FALSE, True/False, 1/0)
strata$FORCED_RACE_REP <- toupper(trimws(as.character(strata$FORCED_RACE_REP))) %in% c("TRUE", "T", "1")

# Controles: si alguno falla, el script se detiene
stopifnot(nrow(geno_all) == 1157, ncol(geno_all) == 20956)
stopifnot(all(unlist(lapply(geno_all, function(x) all(is.na(x) | x %in% c(0, 1, 2))))))
stopifnot(setequal(rownames(geno_all), strata$ID))
stopifnot(sum(strata$STRATUM == "LANDRACE") == 978, sum(strata$STRATUM == "BRED MATERIAL") == 179)
stopifnot(sum(strata$FORCED_RACE_REP) == 22, all(strata$STRATUM[strata$FORCED_RACE_REP] == "LANDRACE"))
cat(sprintf("Datos OK: %d accesiones x %d SNPs, faltantes = %.2f%%\n",
            nrow(geno_all), ncol(geno_all), 100 * mean(is.na(as.matrix(geno_all)))))

# ------------------------------ CORRIDA -------------------------------------
run_stratum <- function(stratum_name, objectives) {
  s     <- strata[strata$STRATUM == stratum_name, ]
  g     <- geno_all[s$ID, , drop = FALSE]
  data  <- coreHunterData(genotypes = genotypes(g, format = "biparental"))
  force <- s$ID[s$FORCED_RACE_REP]
  # Controles previos a Core Hunter: los forzados deben existir en este estrato
  stopifnot(!anyNA(force), all(force %in% rownames(g)))
  stopifnot(length(force) == ifelse(stratum_name == "LANDRACE", 22, 0))
  cat(sprintf("  %s: %d accesiones, %d forzadas, tamano %d\n",
              stratum_name, nrow(g), length(force), SIZE[[stratum_name]]))
  set.seed(SEED)
  core  <- sampleCore(data, obj = objectives, size = SIZE[[stratum_name]],
                      always.selected = if (length(force) > 0) force else integer(0),
                      impr.steps = IMPR_STEPS)
  sel <- core$sel
  stopifnot(length(sel) == SIZE[[stratum_name]], all(force %in% sel), all(sel %in% s$ID))
  data.frame(ID = sel, STRATUM = stratum_name,
             SELECTED_BY = ifelse(sel %in% force, "FORCED_RACE_REP", "CORE HUNTER"),
             stringsAsFactors = FALSE)
}

log <- c(paste("Fecha:", Sys.time()), paste("R:", R.version.string),
         paste("corehunter:", as.character(packageVersion("corehunter"))),
         paste("SEED:", SEED), paste("IMPR_STEPS:", IMPR_STEPS),
         "Exclusion por datos faltantes: ninguna")

for (cfg in names(CONFIGS)) {
  cat("\n=== Configuracion", cfg, "===\n")
  core <- rbind(run_stratum("LANDRACE", CONFIGS[[cfg]]),
                run_stratum("BRED MATERIAL", CONFIGS[[cfg]]))
  stopifnot(nrow(core) == 231, !any(duplicated(core$ID)), sum(core$SELECTED_BY == "FORCED_RACE_REP") == 22)
  write.csv(core, paste0("CORE_STRAT_", cfg, ".csv"), row.names = FALSE)
  cat("Guardado:", paste0("CORE_STRAT_", cfg, ".csv"), "\n")
  log <- c(log, paste("Config", cfg, "OK:", nrow(core), "accesiones"))
}
writeLines(log, "core_estratificada_log.txt")
cat("\nListo. Sube los archivos CORE_STRAT_*.csv y core_estratificada_log.txt para validarlos.\n")
