# =============================================================================
# MetaCiencia: Conocimiento y barreras para la adopción de Reportes Registrados
# en Latinoamérica. Simulación de superpoblación, potencia por Monte Carlo
# y prueba de todos los análisis propuestos.
#
# Versión 1.0.0: esquema de puntajes A (10 ítems de ciencia abierta y 8 de RR;
# numeración P del cuestionario final). Los archivos de results/ se generaron
# con esta versión (doi:10.5281/zenodo.24069486).
# Ejecutar desde la raíz del repositorio: Rscript code/simulacion_analisis_metaciencia.R
# Prueba rápida (3 repeticiones): RRUTA_NSIM=3 Rscript code/simulacion_analisis_metaciencia.R
#
# Estructura
#   0. Paquetes y parámetros
#   1. Generador de datos (superpoblación)
#   2. Depuración (duplicados, respuestas apresuradas, inconsistencias)
#   3. Puntajes e índice de experiencia
#   4. Análisis confirmatorios (C1, C2, C3) con corrección de Holm
#   5. Potencia por Monte Carlo (n = 800, 1200, 1500, 2000)
#   6. Muestra aleatoria de la superpoblación y todos los análisis
#      (descriptivos, exploratorios, DIF, TRI, sensibilidad)
# =============================================================================

# ---- 0. Paquetes y parámetros ------------------------------------------------
library(tidyverse)
library(lme4)       # glmer (C2)
library(glmmTMB)    # beta-binomial mixto (C1)
library(ordinal)    # clmm (C3)
library(mirt)       # TRI 2PL (exploratorio)
library(parallel)
# brms solo se usa si el ajuste frecuentista es singular (ver sección 4)

SEMILLA <- 20261020L
set.seed(SEMILLA)

NSIM       <- as.integer(Sys.getenv("RRUTA_NSIM", "1000"))  # repeticiones por n (RRUTA_NSIM=3 para probar)
OUT        <- "results"      # carpeta de salida
dir.create(OUT, showWarnings = FALSE)
N_GRID     <- c(800, 1200, 1500, 2000)
N_SUPER    <- 20000          # tamaño de la superpoblación simulada
N_MUESTRA  <- 800            # muestra aleatoria para la demostración final
CORES      <- max(1, detectCores() - 1)
ALPHA      <- 0.05           # Holm sobre 3 pruebas confirmatorias
N_CA_ITEMS <- 10             # P21 (9 casillas de contenido) + P23 (1 ítem de selección única)
N_RR_ITEMS <- 8              # P26 (6 casillas de contenido) + P27 + P28
UMBRAL_DESCRIPTIVO_PAIS <- 30
UMBRAL_AREA_OTRA        <- 15

# Efectos verdaderos usados para simular (tamaños "mínimos de interés")
EF <- list(
  beta_exp   = 0.20,   # efecto de la experiencia (z) sobre theta_RR
  gamma_rr   = 0.25,   # efecto del puntaje RR (z) sobre la intención (logit)
  gamma_barr = -0.15,  # efecto del número de barreras (z) sobre la intención
  sd_pais    = 0.30,
  sd_area    = 0.20,
  dif_rr3_pt = -0.8,   # DIF inducido en el ítem RR 3 para la versión en portugués
  a_ca_z     = 0.4,    # efecto de la experiencia (z) sobre theta_CA (vía indirecta hacia RR)
  shift_rr   = 0       # desplazamiento de la familiaridad con RR (se calibra en la sección 4b)
)

paises <- c("Colombia","México","Argentina","Chile","Brasil","Perú",
            "Ecuador","Uruguay","Costa Rica","Bolivia","Venezuela","Otros")
p_paises <- c(.34,.12,.10,.08,.12,.06,.05,.03,.03,.02,.02,.03)
areas <- c("Salud","Psicología","Otras ciencias sociales","Humanidades",
           "Naturales y exactas","Ingeniería y tecnología","Artes y diseño",
           "Agrarias y veterinarias","Educación","Otra")
p_areas <- c(.10,.25,.15,.08,.12,.10,.04,.04,.08,.04)

# ---- 1. Generador de datos ---------------------------------------------------
cut5 <- function(x, cuts) as.integer(cut(x, c(-Inf, cuts, Inf)))

sim_data <- function(n, ef = EF) {
  u_pais <- setNames(rnorm(length(paises), 0, ef$sd_pais), paises)
  u_area <- setNames(rnorm(length(areas), 0, ef$sd_area), areas)

  d <- tibble(
    id      = seq_len(n),
    country = factor(sample(paises, n, TRUE, p_paises), levels = paises),
    area    = factor(sample(areas,  n, TRUE, p_areas),  levels = areas),
    educ    = sample(1:9, n, TRUE, c(.15,.05,.05,.05,.12,.12,.12,.04,.30))
  ) |>
    mutate(
      lang = if_else(country == "Brasil", "pt", "es"),
      has_research = rbinom(n, 1, plogis(0.6 + 0.4 * (educ - 4))),
      z = if_else(has_research == 1, 0.35 * (educ - 5) / 2 + rnorm(n), -1.5),
      # componentes del índice (P7 años, P14 artículos, P9 nivel, P11 liderazgo, P15 roles)
      years  = if_else(has_research == 1, round(pmax(0.5, exp(1.2 + 0.5 * z + rnorm(n, 0, .5)))), 0),
      n_art  = if_else(has_research == 1, rpois(n, exp(0.3 + 0.8 * pmax(z, -1))), 0L),
      nivel_p9   = if_else(has_research == 1, cut5(z + rnorm(n, 0, .7), c(-.5, .6, 1.4))  - 1L, 0L),
      lider_p11   = if_else(has_research == 1, cut5(z + rnorm(n, 0, .9), c(0, 1, 2))       - 1L, 0L),
      n_roles = if_else(has_research == 1, pmin(8L, rpois(n, exp(0.2 + 0.6 * pmax(z, -1)))), 0L),
      # ciencia abierta (latente)
      theta_ca = ef$a_ca_z * z + u_pais[as.character(country)] + u_area[as.character(area)] + rnorm(n),
      fam_ca   = cut5(theta_ca + rlogis(n), c(-1.5, -.5, .5, 1.5)),
      # RR (latente)
      theta_rr = 0.6 * theta_ca + ef$beta_exp * z + rnorm(n, 0, .8),
      fam_rr   = cut5(theta_rr + ef$shift_rr + rlogis(n), c(-.1, .9, 1.8, 3.2)) - 1L,  # 0 = "No"
      rr_answered = fam_rr > 0,
      n_barr  = pmin(8L, rpois(n, exp(1.2 - 0.2 * theta_ca))),
      duration = rlnorm(n, log(7), .4)
    )

  # ítems de ciencia abierta
  b_ca <- seq(-1.2, 1.2, length.out = N_CA_ITEMS)
  for (j in seq_len(N_CA_ITEMS)) {
    d[[sprintf("ca%02d", j)]] <- rbinom(n, 1, plogis(b_ca[j] + d$theta_ca))
  }
  # ítems de RR (solo quienes no respondieron "No" en P25)
  b_rr <- seq(-1, 1, length.out = N_RR_ITEMS)
  for (j in seq_len(N_RR_ITEMS)) {
    shift <- if (j == 3) ifelse(d$lang == "pt", ef$dif_rr3_pt, 0) else 0
    y <- rbinom(n, 1, plogis(b_rr[j] + shift + d$theta_rr))
    d[[sprintf("rr%02d", j)]] <- if_else(d$rr_answered, y, NA_integer_)
  }

  # puntaje observado de RR (cero estructural si "No" en P25), usado para la intención
  rr_obs <- rowSums(d[, sprintf("rr%02d", 1:N_RR_ITEMS)], na.rm = TRUE)
  rr_z   <- as.numeric(scale(rr_obs))
  barr_z <- as.numeric(scale(d$n_barr))
  lat <- ef$gamma_rr * rr_z + ef$gamma_barr * barr_z + 0.2 * d$z +
    u_pais[as.character(d$country)] + u_area[as.character(d$area)] + rlogis(n)
  d$intent_rr <- cut5(lat, c(-1.5, -.5, .5, 1.5))
  d$intent_ca <- cut5(0.8 * lat + 0.3 * d$theta_ca + rlogis(n), c(-1.5, -.5, .5, 1.5))

  d |> select(-theta_ca, -theta_rr, -z)
}

# Contaminación realista: duplicados y respuestas apresuradas al azar
contaminar <- function(d, p_dup = 0.02, p_bot = 0.03) {
  n <- nrow(d)
  bots <- sample(n, round(p_bot * n))
  items <- c(sprintf("ca%02d", 1:N_CA_ITEMS), sprintf("rr%02d", 1:N_RR_ITEMS))
  for (v in items) d[[v]][bots] <- if_else(is.na(d[[v]][bots]), NA_integer_,
                                           rbinom(length(bots), 1, .5))
  d$duration[bots] <- runif(length(bots), 0.5, 1.8)
  dup <- d[sample(n, round(p_dup * n)), ]
  bind_rows(d, dup) |> mutate(id = row_number())
}

# ---- 2. Depuración -----------------------------------------------------------
# Reglas pre-registradas (aplicadas a la base cruda)
depurar <- function(d, mediana_piloto = 7) {
  demo <- c("country","area","educ","lang","has_research","years","n_art","nivel_p9",
            "lider_p11","n_roles")
  items <- c(sprintf("ca%02d", 1:N_CA_ITEMS), sprintf("rr%02d", 1:N_RR_ITEMS),
             "fam_ca","fam_rr","intent_rr","intent_ca","n_barr")
  n0 <- nrow(d)

  # (a) duplicados: mismas respuestas demográficas y >= 90 % de coincidencia en el resto
  d <- d |> group_by(across(all_of(demo))) |> mutate(.g = cur_group_id()) |> ungroup()
  mantener <- rep(TRUE, nrow(d))
  for (g in unique(d$.g[duplicated(d$.g) | duplicated(d$.g, fromLast = TRUE)])) {
    idx <- which(d$.g == g)
    for (k in seq_along(idx)[-1]) {
      for (m in idx[seq_len(k - 1)][mantener[idx[seq_len(k - 1)]]]) {
        a <- unlist(d[m, items]); b <- unlist(d[idx[k], items])
        if (mean(a == b, na.rm = TRUE) >= 0.90) { mantener[idx[k]] <- FALSE; break }
      }
    }
  }
  n_dup <- sum(!mantener)
  d <- d[mantener, ] |> select(-.g)

  # (b) banderas de calidad
  umbral_t <- mediana_piloto / 3
  ca_n <- rowSums(d[, sprintf("ca%02d", 1:N_CA_ITEMS)] == 1, na.rm = TRUE)
  d <- d |> mutate(
    f_tiempo = duration < umbral_t,
    f_logica = (has_research == 0 & (years > 0 | n_art > 0)) |
      (n_art > 15 * pmax(years, 1)),
    f_uniforme = FALSE,     # en datos reales: respuesta idéntica en las escalas (P20, P29, P30)
    f_todas   = ca_n == N_CA_ITEMS,    # proxy; en datos reales: todas las casillas marcadas (incl. "Ninguna") en P21, P24, P26 o P31
    n_debiles = f_uniforme + f_todas,
    excluir = f_tiempo | f_logica | n_debiles >= 2
  )
  resumen <- tibble(inicial = n0, duplicados = n_dup,
                    tiempo = sum(d$f_tiempo), logica = sum(d$f_logica),
                    excluidos_calidad = sum(d$excluir),
                    final = sum(!d$excluir))
  list(datos = filter(d, !excluir), flujo = resumen)
}

# ---- 3. Puntajes e índice de experiencia -------------------------------------
winsor <- function(x, p = 0.99) pmin(x, quantile(x, p, na.rm = TRUE))
re01   <- function(x, mx) x / mx

indice_experiencia <- function(d, version = c("principal", "sin_P15", "pca")) {
  version <- match.arg(version)
  w <- function(x) { x <- winsor(log1p(x)); x / max(x, 1e-9) }
  comp <- tibble(
    anios = w(d$years), arts = w(d$n_art),
    nivel = d$nivel_p9 / 3, lider = d$lider_p11 / 3, roles = pmin(d$n_roles, 8) / 8
  )
  comp[d$has_research == 0, ] <- 0                        # cero estructural
  if (version == "sin_P15") comp <- select(comp, -roles)
  if (version == "pca") {
    pc <- prcomp(comp, scale. = TRUE)$x[, 1]
    if (cor(pc, rowMeans(comp)) < 0) pc <- -pc   # el signo del componente es arbitrario
    return(as.numeric(scale(pc)))
  }
  100 * rowMeans(comp)
}

calcular_puntajes <- function(d) {
  ca <- d[, sprintf("ca%02d", 1:N_CA_ITEMS)]
  rr <- d[, sprintf("rr%02d", 1:N_RR_ITEMS)]
  rr0 <- rr; rr0[is.na(rr0)] <- 0                          # cero por diseño si "No" en P25
  d |> mutate(
    ca_correct = rowSums(ca),
    rr_correct = rowSums(rr0),
    exp_idx    = indice_experiencia(d),
    exp_z      = as.numeric(scale(exp_idx)),
    educ_z     = as.numeric(scale(educ)),
    barr_z     = as.numeric(scale(n_barr)),
    rr_z       = as.numeric(scale(rr_correct)),
    area       = fct_lump_min(area, UMBRAL_AREA_OTRA, other_level = "Otra"),
    intent_f   = factor(intent_rr, levels = 1:5, ordered = TRUE)
  )
}

# ---- 4. Análisis confirmatorios ----------------------------------------------
singular_tmb <- function(m) {
  sds <- unlist(lapply(VarCorr(m)$cond, function(x) attr(x, "stddev")))
  any(sds < 1e-4)
}

# C1: experiencia -> conocimiento de RR (beta-binomial mixto, aleatorios cruzados)
fit_c1 <- function(d) {
  m <- glmmTMB(cbind(rr_correct, N_RR_ITEMS - rr_correct) ~ exp_z + educ_z +
                 (1 | country) + (1 | area), data = d, family = betabinomial)
  cf <- summary(m)$coefficients$cond["exp_z", ]
  list(est = cf[["Estimate"]], p = cf[["Pr(>|z|)"] ], singular = singular_tmb(m), modelo = m)
}

# C2: brecha de familiaridad (CA vs RR), logístico mixto pareado
fit_c2 <- function(d) {
  l <- bind_rows(
    d |> transmute(id, country, area, tipo = "CA", y = as.integer(fam_ca >= 2)),
    d |> transmute(id, country, area, tipo = "RR", y = as.integer(fam_rr >= 1))
  ) |> mutate(tipo = factor(tipo, levels = c("CA", "RR")))
  m <- glmer(y ~ tipo + (1 | id) + (1 | country) + (1 | area), data = l,
             family = binomial, control = glmerControl(optimizer = "bobyqa"))
  cf <- summary(m)$coefficients["tipoRR", ]
  list(est = cf[["Estimate"]], p = cf[["Pr(>|z|)"]], singular = isSingular(m), modelo = m)
}

# C3: conocimiento de RR (y barreras) -> intención de someter RR (clmm)
fit_c3 <- function(d) {
  m <- clmm(intent_f ~ rr_z + barr_z + exp_z + (1 | country) + (1 | area),
            data = d, link = "logit")
  cf <- summary(m)$coefficients
  list(est = cf["rr_z", "Estimate"], p = cf["rr_z", "Pr(>|z|)"],
       est_barr = cf["barr_z", "Estimate"], p_barr = cf["barr_z", "Pr(>|z|)"],
       singular = any(unlist(lapply(ranef(m), function(x) var(x[[1]]))) < 1e-6), modelo = m)
}

confirmatorios <- function(d, completo = FALSE) {
  r1 <- fit_c1(d); r2 <- fit_c2(d); r3 <- fit_c3(d)
  ph <- p.adjust(c(r1$p, r2$p, r3$p), "holm")
  # direcciones pre-registradas: C1 > 0, C2 < 0, C3 > 0
  rej <- c(ph[1] < ALPHA & r1$est > 0, ph[2] < ALPHA & r2$est < 0, ph[3] < ALPHA & r3$est > 0)
  out <- tibble(prueba = c("C1", "C2", "C3"),
                est = c(r1$est, r2$est, r3$est),
                p = c(r1$p, r2$p, r3$p), p_holm = ph, rechaza = rej,
                singular = c(r1$singular, r2$singular, r3$singular))
  if (completo) attr(out, "modelos") <- list(c1 = r1$modelo, c2 = r2$modelo, c3 = r3$modelo)
  out
}

# Respaldo bayesiano si el ajuste es singular (se ejecuta solo si BAYES_RESPALDO)
c1_bayes <- function(d) {
  library(brms)
  brm(rr_correct | trials(N_RR_ITEMS) ~ exp_z + educ_z + (1 | country) + (1 | area),
      data = d |> mutate(N_RR_ITEMS = N_RR_ITEMS), family = beta_binomial(),
      prior = c(prior(normal(0, 1), class = "b"),
                prior(exponential(2), class = "sd")),
      chains = 4, iter = 2000, cores = CORES, seed = 1)
}

# ---- 4b. Calibración de los efectos verdaderos -------------------------------
# Hace que los efectos simulados sean exactamente los efectos mínimos de interés:
#   C1: beta (por 1 DE de experiencia, logit del acierto en ítems RR) = 0,15
#   C2: brecha de familiaridad CA menos RR = 10 puntos porcentuales
#   C3: beta (por 1 DE del puntaje RR, logit acumulado) = 0,15
# Se resuelve con uniroot sobre el promedio de K_CAL muestras (ver abajo).
CALIBRAR <- TRUE
EF_MIN <- list(gap_c2 = 0.10, beta_c1 = 0.15, beta_c3 = 0.15)

# Cada evaluación promedia el efecto estimado sobre K_CAL conjuntos de datos independientes,
# cada uno con sus propios efectos aleatorios de país y área. Una sola muestra, por grande que
# sea, queda calibrada a un sorteo particular de esos efectos (solo hay 12 países), y el efecto
# medio de la potencia (que sortea de nuevo en cada repetición) no coincidía con el objetivo.
# Las semillas son fijas e iguales en cada evaluación, para que uniroot vea una función suave.
N_CAL <- 3000    # tamaño de cada muestra de calibración
K_CAL <- 20      # número de muestras promediadas

calibrar_efectos <- function(n = N_CAL, k = K_CAL, ef = EF, objetivo = EF_MIN) {
  prom <- function(par, valor, f, base) {
    mean(vapply(seq_len(k), function(j) {
      set.seed(base + j)
      d <- sim_data(n, modifyList(ef, setNames(list(valor), par)))
      tryCatch(f(d), error = function(e) NA_real_)
    }, numeric(1)), na.rm = TRUE)
  }
  for (it in 1:2) {   # dos pasadas: los parámetros se afectan entre sí
    gap <- function(s) prom("shift_rr", s,
      function(d) mean(d$fam_ca >= 2) - mean(d$fam_rr >= 1), 1000) - objetivo$gap_c2
    ef$shift_rr <- uniroot(gap, c(-1, 4), tol = 0.002)$root

    est_c1 <- function(b) prom("beta_exp", b,
      function(d) fit_c1(calcular_puntajes(d))$est, 2000) - objetivo$beta_c1
    ef$beta_exp <- uniroot(est_c1, c(-0.8, 0.8), tol = 0.002)$root

    est_c3 <- function(g) prom("gamma_rr", g,
      function(d) fit_c3(calcular_puntajes(d))$est, 3000) - objetivo$beta_c3
    ef$gamma_rr <- uniroot(est_c3, c(0, 1), tol = 0.002)$root
  }
  ef[c("shift_rr", "beta_exp", "gamma_rr")]
}

if (CALIBRAR) {
  # El caché se invalida si cambian el método de calibración, el número de ítems o los efectos mínimos
  cal_file <- file.path(OUT, "calibracion.rds")
  firma <- list("cal_v2", N_CAL, K_CAL, N_CA_ITEMS, N_RR_ITEMS, EF_MIN)
  cache <- if (file.exists(cal_file)) readRDS(cal_file) else NULL
  if (!is.null(cache) && identical(cache$firma, firma)) {
    cal <- cache$cal
  } else {
    cal <- calibrar_efectos(); saveRDS(list(cal = cal, firma = firma), cal_file)
  }
  EF <- modifyList(EF, cal)
  cat("\n== Efectos calibrados ==\n"); print(unlist(cal))
}

# ---- 5. Potencia por Monte Carlo ---------------------------------------------
una_rep <- function(n, ef = EF) {
  tryCatch({
    d <- sim_data(n, ef) |> calcular_puntajes()
    confirmatorios(d) |> mutate(n = n, error = NA_character_)
  }, error = function(e) tibble(prueba = c("C1","C2","C3"), est = NA_real_, p = NA_real_,
                                p_holm = NA_real_, rechaza = NA, singular = NA, n = n,
                                error = conditionMessage(e)))   # el mensaje queda en el crudo
}

potencia <- function(ns = N_GRID, nsim = NSIM, ef = EF) {
  # Se conservan los resultados de cada repetición (p, estimación, rechazo) en
  # attr(, "crudo"), para graficar su distribución sin volver a simular.
  crudo <- map_dfr(ns, function(n) {
    # semilla propia por repetición: resultado idéntico sin importar el número de núcleos
    mclapply(seq_len(nsim), function(i) {
      set.seed(SEMILLA + as.integer(n) * 1000L + i)
      una_rep(n, ef)
    }, mc.cores = CORES) |>
      bind_rows(.id = "rep")
  })
  res <- crudo |> group_by(n, prueba) |>
    summarise(reps_validas = sum(!is.na(rechaza)),
              potencia = mean(rechaza, na.rm = TRUE),
              est_media = mean(est, na.rm = TRUE),
              prop_singular = mean(singular, na.rm = TRUE), .groups = "drop")
  attr(res, "crudo") <- crudo
  res
}

EF_NULO <- modifyList(EF, list(beta_exp = 0, a_ca_z = 0, gamma_rr = 0, gamma_barr = 0))

cat("\n== Potencia (Monte Carlo) ==\n")
res_potencia <- potencia()
print(res_potencia, n = Inf)
write_csv(res_potencia, file.path(OUT, "resultados_potencia.csv"))

# ---- 5a. Error tipo I: con efectos nulos en C1 y C3 (C2 no es nulo por construcción) ----
# Tasa de rechazo esperada: no mayor que 0,025 (solo cuentan los rechazos en la dirección
# hipotetizada); con C2 rechazada primero, Holm exige menos a C1 y C3, así que puede ser menor.
NSIM_NULA <- as.integer(Sys.getenv("RRUTA_NSIM_NULA", "500"))
cat("\n== Error tipo I (efectos nulos, n = 1500) ==\n")
res_nula <- potencia(ns = 1500, nsim = NSIM_NULA, ef = EF_NULO) |>
  filter(prueba != "C2") |> rename(tasa_rechazo = potencia)
print(res_nula)
write_csv(res_nula, file.path(OUT, "resultados_error_tipo1.csv"))

# ---- 5b. Figura: distribución de los valores p de la simulación de potencia ----
# Usa los resultados ya simulados (no agrega cálculos). Cada curva es la función de
# distribución acumulada de los valores p (ajustados por Holm) de las repeticiones.
# Su altura en la línea vertical (α = 0,05) es la potencia: cuanto más a la izquierda
# y arriba esté la curva, más potencia. Un valor p cuenta solo si el efecto estimado
# va en la dirección pre-registrada (C1 > 0, C2 < 0, C3 > 0); si no, se cuenta como 1.
crudo <- attr(res_potencia, "crudo")
saveRDS(crudo, file.path(OUT, "resultados_potencia_crudo.rds"))

dir_esp  <- c(C1 = 1, C2 = -1, C3 = 1)
nom_c    <- c(C1 = "C1. Experiencia y conocimiento de RR",
              C2 = "C2. Brecha de familiaridad (CA frente a RR)",
              C3 = "C3. Conocimiento de RR e intención")
niveles_n <- sort(unique(crudo$n))

pdat <- crudo |>
  filter(!is.na(p_holm)) |>
  mutate(p_dir = if_else(sign(est) == unname(dir_esp[prueba]), p_holm, 1),
         p_dir = pmax(p_dir, 1e-10),                       # evita log(0); se acumula en el borde
         n_f = factor(n, levels = niveles_n, labels = format(niveles_n, big.mark = ".", decimal.mark = ",", trim = TRUE)),
         prueba_f = factor(prueba, levels = names(nom_c), labels = nom_c))

pot_pts <- pdat |> group_by(prueba_f, n_f) |>
  summarise(pot = mean(p_dir < ALPHA), .groups = "drop")

p5 <- ggplot(pdat, aes(p_dir, color = n_f)) +
  stat_ecdf(geom = "step", linewidth = 1, pad = FALSE) +
  geom_vline(xintercept = ALPHA, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "grey40") +
  geom_point(data = pot_pts, aes(ALPHA, pot, color = n_f), size = 2.4, show.legend = FALSE) +
  facet_wrap(~ prueba_f, nrow = 1) +
  scale_x_log10(breaks = 10^seq(-10, 0, 2), labels = scales::label_log(), limits = c(1e-10, 1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_viridis_d(end = 0.85, name = "Tamaño de muestra (n)") +
  labs(x = "Valor p ajustado por Holm (escala logarítmica)",
       y = "Proporción de simulaciones con p ≤ x",
       title = "Distribución de los valores p en la simulación de potencia",
       subtitle = "La altura de cada curva en la línea punteada vertical (α = 0,05) es la potencia; la línea horizontal marca 80 %",
       caption = paste0(max(crudo$rep |> as.integer()), " repeticiones por tamaño de muestra. ",
                        "Los valores p menores a 10⁻¹⁰ se acumulan en el borde izquierdo.")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 10))
ggsave(file.path(OUT, "fig_potencia.png"), p5, width = 14, height = 5.5, dpi = 300, bg = "white")

# ---- 6. Superpoblación, muestra aleatoria y análisis completo ----------------
set.seed(SEMILLA + 1L)
super <- sim_data(N_SUPER) |> contaminar()
muestra <- super |> slice_sample(n = N_MUESTRA)

dep <- depurar(muestra)
cat("\n== Flujo de participantes ==\n"); print(dep$flujo)
datos <- calcular_puntajes(dep$datos)

# 6.1 Descriptivos: proporciones con IC de Wilson, por país (n >= umbral)
wilson <- function(x, n) {
  ic <- prop.test(x, n, correct = FALSE)$conf.int
  tibble(prop = x / n, lwr = ic[1], upr = ic[2])
}
desc_pais <- datos |>
  count(country, name = "n") |>
  left_join(datos |> group_by(country) |>
              summarise(x_rr = sum(fam_rr >= 1), x_ca = sum(fam_ca >= 2)), by = "country") |>
  filter(n >= UMBRAL_DESCRIPTIVO_PAIS) |>
  rowwise() |>
  mutate(rr = list(wilson(x_rr, n)), ca = list(wilson(x_ca, n))) |>
  unnest_wider(rr, names_sep = "_") |> unnest_wider(ca, names_sep = "_")
cat("\n== Familiaridad por país (n >= 30), IC Wilson 95% ==\n"); print(desc_pais)

# 6.2 Confiabilidad (KR-20) de los puntajes
kr20 <- function(m) {
  m <- as.matrix(m); m[is.na(m)] <- 0
  k <- ncol(m); p <- colMeans(m); v <- var(rowSums(m))
  (k / (k - 1)) * (1 - sum(p * (1 - p)) / v)
}
cat("\nKR-20 CA:", round(kr20(datos[, sprintf("ca%02d", 1:N_CA_ITEMS)]), 3),
    " | KR-20 RR:", round(kr20(datos[, sprintf("rr%02d", 1:N_RR_ITEMS)]), 3), "\n")

# 6.3 Índice de experiencia: intercorrelaciones y VIF
comp <- datos |> transmute(anios = log1p(years), arts = log1p(n_art), nivel = nivel_p9,
                           lider = lider_p11, roles = n_roles)
cat("\n== Correlaciones entre componentes del índice ==\n"); print(round(cor(comp), 2))
vif_simple <- function(df) sapply(names(df), function(v)
  1 / (1 - summary(lm(reformulate(setdiff(names(df), v), v), data = df))$r.squared))
cat("\nVIF exp_z, educ_z:\n"); print(round(vif_simple(datos |> select(exp_z, educ_z)), 2))

# 6.4 Confirmatorios con los datos de la muestra
conf <- confirmatorios(datos, completo = TRUE)
cat("\n== Pruebas confirmatorias (Holm) ==\n"); print(conf)

# 6.4b Figuras de la muestra: un panel por prueba confirmatoria
# (las curvas son las predicciones de los modelos con efectos aleatorios en cero;
#  los puntos son los datos observados de la muestra)
pal <- c(CA = "#0072B2", RR = "#D55E00")            # paleta apta para daltonismo
col_rr <- unname(pal["RR"])
fmt_p <- function(p) if (p < 0.001) "p < 0,001" else paste0("p = ", format(round(p, 3), decimal.mark = ","))
sub_c <- function(k) {
  f <- conf[conf$prueba == k, ]
  paste0("β = ", format(round(f$est, 2), decimal.mark = ","), "; ", fmt_p(f$p_holm), " (Holm)")
}
tema <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 10))

# Panel C1: experiencia -> aciertos en ítems de RR
m1 <- attr(conf, "modelos")$c1
mu_e <- mean(datos$exp_idx); sd_e <- sd(datos$exp_idx)
g1 <- tibble(exp_idx = seq(min(datos$exp_idx), max(datos$exp_idx), length.out = 100)) |>
  mutate(exp_z = (exp_idx - mu_e) / sd_e, educ_z = 0,
         country = datos$country[1], area = datos$area[1])
pr1 <- predict(m1, newdata = g1, re.form = NA, se.fit = TRUE, type = "link",
               allow.new.levels = TRUE)
g1 <- g1 |> mutate(fit = plogis(pr1$fit),
                   lwr = plogis(pr1$fit - 1.96 * pr1$se.fit),
                   upr = plogis(pr1$fit + 1.96 * pr1$se.fit))
p1 <- ggplot() +
  geom_jitter(data = datos, aes(exp_idx, rr_correct / N_RR_ITEMS),
              width = 0, height = 0.015, alpha = 0.25, size = 1, color = "grey40") +
  geom_ribbon(data = g1, aes(exp_idx, ymin = lwr, ymax = upr), fill = col_rr, alpha = 0.2) +
  geom_line(data = g1, aes(exp_idx, fit), color = col_rr, linewidth = 1.1) +
  scale_y_continuous(labels = scales::percent) + coord_cartesian(ylim = c(0, 1)) +
  labs(title = "C1. Experiencia y conocimiento de RR", subtitle = sub_c("C1"),
       x = "Índice de experiencia investigativa (0 a 100)",
       y = "Aciertos en los ítems de RR") + tema

# Panel C2: familiaridad con CA frente a RR (total y por país)
tot <- bind_rows(
  bind_cols(tibble(grupo = "Total", tipo = "CA"), wilson(sum(datos$fam_ca >= 2), nrow(datos))),
  bind_cols(tibble(grupo = "Total", tipo = "RR"), wilson(sum(datos$fam_rr >= 1), nrow(datos)))
)
por_pais <- desc_pais |>
  ungroup() |>
  transmute(grupo = as.character(country), ca_prop, ca_lwr, ca_upr, rr_prop, rr_lwr, rr_upr) |>
  pivot_longer(-grupo, names_to = c("tipo", ".value"), names_sep = "_") |>
  mutate(tipo = toupper(tipo))
dat2 <- bind_rows(tot, por_pais) |>
  mutate(grupo = factor(grupo, levels = c("Total", as.character(desc_pais$country))))
p2 <- ggplot(dat2, aes(prop, fct_rev(grupo), color = tipo)) +
  geom_pointrange(aes(xmin = lwr, xmax = upr), position = position_dodge(width = 0.55), size = 0.4) +
  scale_color_manual(values = pal, name = NULL,
                     labels = c(CA = "Ciencia abierta (P20 ≥ 2)", RR = "Reportes Registrados (P25 ≠ No)")) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "C2. Familiaridad: ciencia abierta frente a RR", subtitle = sub_c("C2"),
       x = "Proporción con familiaridad (IC de Wilson 95 %)", y = NULL) +
  tema + theme(legend.position = "bottom", legend.direction = "vertical")

# Panel C3: conocimiento de RR -> probabilidad de intención alta (P30 = 4 o 5)
m3 <- attr(conf, "modelos")$c3
cf3 <- coef(m3)
mu_r <- mean(datos$rr_correct); sd_r <- sd(datos$rr_correct)
g3 <- tibble(rr = seq(0, N_RR_ITEMS, length.out = 100)) |>
  mutate(rr_z = (rr - mu_r) / sd_r,
         fit = 1 - plogis(cf3[["3|4"]] - cf3[["rr_z"]] * rr_z),
         lwr = NA_real_, upr = NA_real_)
vc3 <- tryCatch(vcov(m3), error = function(e) NULL)
if (!is.null(vc3) && all(c("3|4", "rr_z") %in% rownames(vc3))) {
  set.seed(1)
  dr <- MASS::mvrnorm(2000, mu = cf3[c("3|4", "rr_z")], Sigma = vc3[c("3|4", "rr_z"), c("3|4", "rr_z")])
  bandas <- t(sapply(g3$rr_z, function(z) quantile(1 - plogis(dr[, 1] - dr[, 2] * z), c(0.025, 0.975))))
  g3$lwr <- bandas[, 1]; g3$upr <- bandas[, 2]
}
obs3 <- datos |> group_by(rr_correct) |> summarise(n = n(), p = mean(intent_rr >= 4), .groups = "drop")
p3 <- ggplot() +
  geom_ribbon(data = g3, aes(rr, ymin = lwr, ymax = upr), fill = col_rr, alpha = 0.2) +
  geom_line(data = g3, aes(rr, fit), color = col_rr, linewidth = 1.1) +
  geom_point(data = obs3, aes(rr_correct, p, size = n), color = "grey30", alpha = 0.7) +
  scale_size_area(max_size = 5, name = "Personas") +
  scale_y_continuous(labels = scales::percent) + coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, N_RR_ITEMS, 2)) +
  labs(title = "C3. Conocimiento de RR e intención de someter un RR", subtitle = sub_c("C3"),
       x = paste0("Puntaje de conocimiento de RR (0 a ", N_RR_ITEMS, ")"),
       y = "Intención alta (P30 = 4 o 5)") + tema + theme(legend.position = "bottom")

ggsave(file.path(OUT, "fig_C1.png"), p1, width = 6, height = 5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig_C2.png"), p2, width = 6.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig_C3.png"), p3, width = 6, height = 5, dpi = 300, bg = "white")
if (requireNamespace("patchwork", quietly = TRUE)) {
  fig <- patchwork::wrap_plots(p1, p2, p3, nrow = 1) +
    patchwork::plot_annotation(
      tag_levels = "A",
      caption = paste0("Datos simulados: muestra aleatoria de n = ", nrow(datos),
                       " tras la depuración. Curvas: predicciones de los modelos con efectos aleatorios de país y área en cero."))
  ggsave(file.path(OUT, "fig_C1_C2_C3.png"), fig, width = 18, height = 6.5, dpi = 300, bg = "white")
} else message("Instala 'patchwork' para obtener la figura combinada; se guardaron las tres figuras por separado.")

# 6.5 Sensibilidad
# (a) índice alternativo (PCA y sin P15)
for (v in c("pca", "sin_P15")) {
  dd <- datos |> mutate(exp_z = as.numeric(scale(indice_experiencia(datos, v))))
  cat("\nSensibilidad C1 con índice", v, ":\n")
  print(round(summary(fit_c1(dd)$modelo)$coefficients$cond["exp_z", ], 4))
}
# (b) solo quienes contestaron el bloque de RR (P25 != "No")
cat("\nSensibilidad C1 y C3 solo con P25 != 'No':\n")
print(confirmatorios(filter(datos, rr_answered)) |> filter(prueba != "C2"))
# (c) área como efecto fijo
m_fix <- glmmTMB(cbind(rr_correct, N_RR_ITEMS - rr_correct) ~ exp_z + educ_z + area +
                   (1 | country), data = datos, family = betabinomial)
cat("\nC1 con área como efecto fijo:\n"); print(round(summary(m_fix)$coefficients$cond["exp_z", ], 4))
# (d) con y sin depuración
conf_cruda <- confirmatorios(calcular_puntajes(muestra))
cat("\nC1-C3 sin depuración:\n"); print(conf_cruda)

# 6.6 Exploratorio: ICC (varianza por país y área) del modelo de C1
vc <- VarCorr(attr(conf, "modelos")$c1)$cond
cat("\nDesviaciones estándar aleatorias C1:\n"); print(lapply(vc, attr, "stddev"))

# 6.7 DIF por versión de idioma (regresión logística, Swaminathan y Rogers)
nagelkerke_delta <- function(m0, m1) {
  n <- nobs(m1); r2 <- function(m) (1 - exp((m$deviance - m$null.deviance) / n)) /
    (1 - exp(-m$null.deviance / n))
  r2(m1) - r2(m0)
}
dif_res <- map_dfr(sprintf("rr%02d", 1:N_RR_ITEMS), function(it) {
  dd <- datos |> filter(rr_answered, !is.na(.data[[it]])) |>
    mutate(resto = rr_correct - .data[[it]], pt = as.integer(lang == "pt"))
  if (sum(dd$pt) < 20) return(tibble(item = it, delta_R2 = NA_real_, clase = "sin datos"))
  y <- dd[[it]]
  m0 <- glm(y ~ resto, data = dd, family = binomial)
  m1 <- glm(y ~ resto * pt, data = dd, family = binomial)
  d2 <- nagelkerke_delta(m0, m1)
  tibble(item = it, delta_R2 = d2,
         clase = case_when(d2 < 0.035 ~ "despreciable", d2 < 0.07 ~ "moderado", TRUE ~ "grande"))
})
cat("\n== DIF (español vs portugués), ítems RR; ítem 3 tiene DIF inducido ==\n")
print(dif_res)

# 6.8 TRI 2PL conjunta (exploratorio): dificultades de ítems CA vs RR
items_irt <- datos |> filter(rr_answered) |>
  select(matches("^(ca|rr)[0-9]{2}$")) |> select(where(~ var(.x, na.rm = TRUE) > 0))
irt <- tryCatch(mirt(items_irt, 1, itemtype = "2PL", verbose = FALSE), error = function(e) NULL)
if (!is.null(irt)) {
  par <- coef(irt, IRTpars = TRUE, simplify = TRUE)$items
  cat("\nTRI 2PL: dificultad media CA =", round(mean(par[grepl("^ca", rownames(par)), "b"]), 2),
      " | RR =", round(mean(par[grepl("^rr", rownames(par)), "b"]), 2), "\n")
}

# 6.9 Descriptivos de barreras y razones (ejemplo: número de barreras por nivel)
cat("\nBarreras medias por nivel de formación:\n")
print(datos |> group_by(educ) |> summarise(n = n(), barreras = mean(n_barr)) |> arrange(educ))

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("\nListo. Resultados en results/\n")
