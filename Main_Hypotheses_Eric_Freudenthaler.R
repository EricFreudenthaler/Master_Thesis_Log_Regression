library(ordinal)
library(ggplot2)
library(ggeffects)
library(marginaleffects)
library(modelsummary)

daten <- read.csv("/Users/ericfreudenthaler/Documents/MasterArbeit/Projects/R Regression/CSV_Master_Thesis_V2.csv", #wird noch ausgetauscht falls ich es auf Github laden soll
                  stringsAsFactors = FALSE, #keeps text as text 
                  na.strings = c("", "NA", "n.a.", "#N/A", "#DIV/0!", "#VALUE!")) #these values should be treated as missing data

#Dependent Variable is set in the correct order from high (Withdrawal) to low (Digging In)
daten$CELI.Grade <- factor(
  daten$CELI.Grade,
  levels = c("Digging In", "Buying Time", "Scaling Back", "Suspension", "Withdrawal"),
  ordered = TRUE
)

#Preparing Controls into the right values, and adds +1 to prevent log(0), the log compresses thes scale and adresses the right skewedness
prepare_controls <- function(df) {
  df$log_firm_size <- log(as.numeric(df$Global_Revenue_2021_m) + 1)
  df$ROA <- as.numeric(df$ROA.using.Profit...Loss..before.tax.2021....)
  df$industry_grouped <- df$GICS.Industry.Sector
  df$industry_grouped[df$industry_grouped %in% c("Real Estate", "Utilities")] <- "Other"
  df$industry_grouped <- factor(df$industry_grouped)
  df$industry_grouped <- relevel(df$industry_grouped, ref = "Industrials")
  df$region <- ifelse(df$Country.CELI == "United States", "US",
               ifelse(df$Country.CELI == "United Kingdom", "UK",
               ifelse(df$Country.CELI == "Germany", "Germany",
               ifelse(df$Country.CELI == "France", "France",
               ifelse(df$Country.CELI == "Finland", "Finland", "Other EU")))))
  df$region <- factor(df$region)
  df$region <- relevel(df$region, ref = "US")
  return(df) #Returns dataframe with all four control variables
}

#Diagnostics
run_diagnostics <- function(modell_full, df, label) {
  cat("\n=== Diagnostics:", label, "===\n")
  # POM Test takes the Proportional Odds Model for every predictor: if p < 0.05, that predictor violates the POM
  cat("\nNominal Test (POM):\n")
  print(nominal_test(modell_full))
  # Sparsity: Checks the threshold parameters when n is very low as this can affect the models stability
  cat("\nCELI Verteilung", label, ":\n")
  print(table(df$CELI.Grade))
  print(round(prop.table(table(df$CELI.Grade)) * 100, 1))
  # LR Test: This creates the null model: no predictors, only the thresholds (the 1 gives every firm has the same probability distribution)
  modell_null <- clm(CELI.Grade ~ 1, data = modell_full$model)
  cat("\nLR Test:\n")
  #Compares the log-likelihood of the null model vs. the full model (if p<0.05 model with predictors is sign. better than null mo)
  print(anova(modell_null, modell_full))
  # Nagelkerke Pseudo R2: First the loglogik extracts the likelihood value (higher=better fit)
  #Nagelkerke model measures how much of the maximum % f the maximum possible improvement over the null model can be captured
  ll_full <- as.numeric(logLik(modell_full))
  ll_null <- as.numeric(logLik(modell_null))
  n <- nobs(modell_full)
  nagelkerke_r2 <- (1 - exp((2/n) * (ll_null - ll_full))) / (1 - exp((2/n) * ll_null))
  cat("Nagelkerke Pseudo R²:", round(nagelkerke_r2, 4), "\n") #gerundet auf 4 komma stellen
  # AIC and BIC of the full model (lower = better; BIC penalises extra parameters more strongly)
  cat("AIC:", round(AIC(modell_full), 2), " BIC:", round(BIC(modell_full), 2), "\n")
}

# Computes 95% confidence intervals for every coefficient with profile likelihood. 
odds_ratios <- function(modell) {
  ci <- confint(modell)
  or_table <- data.frame(
    Variable = rownames(ci),
    OR = round(exp(coef(modell)[rownames(ci)]), 4),
    CI_lower = round(exp(ci[, 1]), 4),
    CI_upper = round(exp(ci[, 2]), 4)
  )
  return(or_table)
}

# Regression H1a: Revenue Share (tbd in Results Chapter 6)
daten$Revenue_Share_RU_2021_ <- as.numeric(daten$Revenue_Share_RU_2021_)
daten_rev <- daten[!is.na(daten$Revenue_Share_RU_2021_) & daten$Revenue_Share_RU_2021_ <= 1, ] #Filter to remove revenue share abover 100%
cat("H1a Beobachtungen:", nrow(daten_rev), "\n")
print(summary(daten_rev$Revenue_Share_RU_2021_))

daten_rev <- prepare_controls(daten_rev)
daten_rev$rev_share_pct <- daten_rev$Revenue_Share_RU_2021_ * 100 #Revenue shares rescaled to percentage points

# Base Model shows the raw effect before adding controls
modell_base_h1a <- clm(CELI.Grade ~ rev_share_pct, data = daten_rev)
summary(modell_base_h1a)

#Full Modell
modell_h1a <- clm(CELI.Grade ~ rev_share_pct + log_firm_size + industry_grouped + ROA + region,
                  data = daten_rev)
summary(modell_h1a)

run_diagnostics(modell_h1a, daten_rev, "H1a")
print(odds_ratios(modell_h1a))
ame_h1a <- avg_slopes(modell_h1a, variables = "rev_share_pct")
print(ame_h1a)


#H1b: Asset share (tbd in Results Chapter 6)
daten$Asset_Share_RU_2021 <- as.numeric(daten$Asset_Share_RU_2021)
daten_asset <- daten[!is.na(daten$Asset_Share_RU_2021) & daten$Asset_Share_RU_2021 <= 1, ]
cat("H1b Beobachtungen:", nrow(daten_asset), "\n")
print(summary(daten_asset$Asset_Share_RU_2021))

daten_asset <- prepare_controls(daten_asset)
daten_asset$asset_share_pct <- daten_asset$Asset_Share_RU_2021 * 100

#Base Model shows the raw effect before adding controls
modell_base_h1b <- clm(CELI.Grade ~ asset_share_pct, data = daten_asset)
summary(modell_base_h1b)

#Full Modell
modell_h1b <- clm(CELI.Grade ~ asset_share_pct + log_firm_size + industry_grouped + ROA + region,
                  data = daten_asset)
summary(modell_h1b)

run_diagnostics(modell_h1b, daten_asset, "H1b")
print(odds_ratios(modell_h1b))
ame_h1b <- avg_slopes(modell_h1b, variables = "asset_share_pct")
print(ame_h1b)

# H1c: Employee Share of Worldwide Workforce (tbd in Results Chapter 6)
daten$Employee_Share_in_Russia <- as.numeric(daten$Employee_Share_in_Russia)
daten_emp <- daten[!is.na(daten$Employee_Share_in_Russia) & daten$Employee_Share_in_Russia <= 1, ]
cat("H1c Beobachtungen:", nrow(daten_emp), "\n")
print(summary(daten_emp$Employee_Share_in_Russia))

daten_emp <- prepare_controls(daten_emp)
daten_emp$emp_share_pct <- daten_emp$Employee_Share_in_Russia * 100

#Base Model
modell_base_h1c <- clm(CELI.Grade ~ emp_share_pct, data = daten_emp)
summary(modell_base_h1c)

#Full Model
modell_h1c <- clm(CELI.Grade ~ emp_share_pct + log_firm_size + industry_grouped + ROA + region,
                  data = daten_emp)
summary(modell_h1c)

run_diagnostics(modell_h1c, daten_emp, "H1c")
print(odds_ratios(modell_h1c))
ame_h1c <- avg_slopes(modell_h1c, variables = "emp_share_pct")
print(ame_h1c)

# combined standardized index is non-significant


z <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)

daten$rev_sh <- as.numeric(daten$Revenue_Share_RU_2021_)
daten$ass_sh <- as.numeric(daten$Asset_Share_RU_2021)
daten$emp_sh <- as.numeric(daten$Employee_Share_in_Russia)

# firms with all three shares valid (and <= 1) and a CELI grade
daten_exp <- daten[!is.na(daten$rev_sh) & daten$rev_sh <= 1 &
                     !is.na(daten$ass_sh) & daten$ass_sh <= 1 &
                     !is.na(daten$emp_sh) & daten$emp_sh <= 1 &
                     !is.na(daten$CELI.Grade), ]

daten_exp$exp_index <- rowMeans(cbind(z(daten_exp$rev_sh),
                                      z(daten_exp$ass_sh),
                                      z(daten_exp$emp_sh)))

daten_exp <- prepare_controls(daten_exp)

modell_exp <- clm(CELI.Grade ~ exp_index + log_firm_size + industry_grouped + ROA + region,
                  data = daten_exp)
summary(modell_exp)
run_diagnostics(modell_exp, daten_exp, "Composite economic exposure")
print(odds_ratios(modell_exp))
print(avg_slopes(modell_exp, variables = "exp_index"))


#H2: Sanction Rating (tbd in Results Chapter 6) 
daten$Sanction_Rating <- as.numeric(daten$Sanction_Rating)
daten_sanc <- daten[!is.na(daten$Sanction_Rating) & !is.na(daten$CELI.Grade), ]
daten_sanc$sanction <- factor(daten_sanc$Sanction_Rating) # treats 0, 1, 2 as categories
cat("H2 Beobachtungen:", nrow(daten_sanc), "\n")
cat("Sanction Verteilung:\n")
print(table(daten_sanc$sanction))

# Controls without industry (because industry is nearly collinear with sanction rating)
daten_sanc$log_firm_size <- log(as.numeric(daten_sanc$Global_Revenue_2021_m) + 1)
daten_sanc$ROA <- as.numeric(daten_sanc$ROA.using.Profit...Loss..before.tax.2021....)
daten_sanc$region <- ifelse(daten_sanc$Country.CELI == "United States", "US",
                            ifelse(daten_sanc$Country.CELI == "United Kingdom", "UK",
                                   ifelse(daten_sanc$Country.CELI == "Germany", "Germany",
                                          ifelse(daten_sanc$Country.CELI == "France", "France",
                                                 ifelse(daten_sanc$Country.CELI == "Finland", "Finland", "Other EU")))))
daten_sanc$region <- factor(daten_sanc$region)
daten_sanc$region <- relevel(daten_sanc$region, ref = "US")

# Basismodell
modell_base_h2 <- clm(CELI.Grade ~ sanction, data = daten_sanc)
summary(modell_base_h2)

# Volles Modell (NO industry control)
modell_h2 <- clm(CELI.Grade ~ sanction + log_firm_size + ROA + region, #no industry
                 data = daten_sanc)
summary(modell_h2)

run_diagnostics(modell_h2, daten_sanc, "H2")
print(odds_ratios(modell_h2))
ame_h2 <- avg_slopes(modell_h2, variables = "sanction")
print(ame_h2)

# H3a: Subsidiaries in Russia (tbd in Results Chapter 6)
daten$Subsidiaries_RU_Latest <- as.numeric(daten$Subsidiaries_RU_Latest)
daten_sub <- daten[!is.na(daten$Subsidiaries_RU_Latest) & !is.na(daten$CELI.Grade), ]
cat("H3a Beobachtungen:", nrow(daten_sub), "\n")
print(summary(daten_sub$Subsidiaries_RU_Latest))

daten_sub$log_subs <- log(daten_sub$Subsidiaries_RU_Latest + 1) #log transformation because of righ skewedness of distribution
daten_sub <- prepare_controls(daten_sub)

#Base Model
modell_base_h3a <- clm(CELI.Grade ~ log_subs, data = daten_sub)
summary(modell_base_h3a)

#Full Model
modell_h3a <- clm(CELI.Grade ~ log_subs + log_firm_size + industry_grouped + ROA + region,
                  data = daten_sub)
summary(modell_h3a)

run_diagnostics(modell_h3a, daten_sub, "H3a")
print(odds_ratios(modell_h3a))
ame_h3a <- avg_slopes(modell_h3a, variables = "log_subs")
print(ame_h3a)

# H3b: Employees in Russian subsidiaries: Industry controls excluded due to complete separation because some industry × CELI combinations have zero firms.

daten_emp_ru <- daten[!is.na(as.numeric(daten$Employees_Subsidiaries_RU)), ]
daten_emp_ru$emp_ru <- as.numeric(daten_emp_ru$Employees_Subsidiaries_RU)
cat("H3b Beobachtungen:", nrow(daten_emp_ru), "von", nrow(daten), "\n")
cat("Coverage:", round(nrow(daten_emp_ru)/nrow(daten)*100, 1), "%\n\n")
print(summary(daten_emp_ru$emp_ru))

daten_emp_ru$log_emp_ru <- log(daten_emp_ru$emp_ru + 1)
daten_emp_ru$log_firm_size <- log(as.numeric(daten_emp_ru$Global_Revenue_2021_m) + 1)
daten_emp_ru$ROA <- as.numeric(daten_emp_ru$ROA.using.Profit...Loss..before.tax.2021....)
daten_emp_ru$region <- ifelse(daten_emp_ru$Country.CELI == "United States", "US",
                              ifelse(daten_emp_ru$Country.CELI == "United Kingdom", "UK",
                                     ifelse(daten_emp_ru$Country.CELI == "Germany", "Germany",
                                            ifelse(daten_emp_ru$Country.CELI == "France", "France",
                                                   ifelse(daten_emp_ru$Country.CELI == "Finland", "Finland", "Other EU")))))
daten_emp_ru$region <- factor(daten_emp_ru$region)
daten_emp_ru$region <- relevel(daten_emp_ru$region, ref = "US")

#Base Model
modell_base_h3b <- clm(CELI.Grade ~ log_emp_ru, data = daten_emp_ru)
summary(modell_base_h3b)

#Primary Model (without industry due to complete separation)
modell_h3b <- clm(CELI.Grade ~ log_emp_ru + log_firm_size + ROA + region, data = daten_emp_ru)
summary(modell_h3b)

run_diagnostics(modell_h3b, daten_emp_ru, "H3b")
print(odds_ratios(modell_h3b))
ame_h3b <- avg_slopes(modell_h3b, variables = "log_emp_ru")
print(ame_h3b)

# H3c: Assets in Russian subsidiaries (tbd in Results Chapter 6)
daten$Assets_RU_Subsidiaries <- as.numeric(daten$Assets_RU_Subsidiaries)
daten_assets_ru <- daten[!is.na(daten$Assets_RU_Subsidiaries) & !is.na(daten$CELI.Grade), ]
cat("H3c Beobachtungen:", nrow(daten_assets_ru), "von", nrow(daten), "\n")
cat("Coverage:", round(nrow(daten_assets_ru)/nrow(daten)*100, 1), "%\n\n")
print(summary(daten_assets_ru$Assets_RU_Subsidiaries))

daten_assets_ru$log_assets_ru <- log(daten_assets_ru$Assets_RU_Subsidiaries + 1) #log transformation because of right skewedness
daten_assets_ru <- prepare_controls(daten_assets_ru)

#Base Model
modell_base_h3c <- clm(CELI.Grade ~ log_assets_ru, data = daten_assets_ru)
summary(modell_base_h3c)

#Primary Model (without industry, matching H3b: dropped due to complete/near-complete
#separation in the small subsample, i.e. empty industry x CELI cells)
modell_h3c <- clm(CELI.Grade ~ log_assets_ru + log_firm_size + ROA + region,
                  data = daten_assets_ru)
summary(modell_h3c)

run_diagnostics(modell_h3c, daten_assets_ru, "H3c")
print(odds_ratios(modell_h3c))
ame_h3c <- avg_slopes(modell_h3c, variables = "log_assets_ru")
print(ame_h3c)

# H3c robustness check: model WITH industry control (reported for transparency)
modell_h3c_ind <- clm(CELI.Grade ~ log_assets_ru + log_firm_size + industry_grouped + ROA + region,
                      data = daten_assets_ru)
summary(modell_h3c_ind)
run_diagnostics(modell_h3c_ind, daten_assets_ru, "H3c with industry")
print(odds_ratios(modell_h3c_ind))
ame_h3c_ind <- avg_slopes(modell_h3c_ind, variables = "log_assets_ru")
print(ame_h3c_ind)

# H4a: Cash/Asset Ratio (tbd in Results Chapter 6)
daten_cash <- daten[!is.na(as.numeric(daten$CashAsset.Ratio_2021)), ]
daten_cash$cash_asset <- as.numeric(daten_cash$CashAsset.Ratio_2021)
cat("H4a Beobachtungen:", nrow(daten_cash), "von", nrow(daten), "\n")
cat("Coverage:", round(nrow(daten_cash)/nrow(daten)*100, 1), "%\n\n")
print(summary(daten_cash$cash_asset))

daten_cash <- prepare_controls(daten_cash)

#Base Model
modell_base_h4a <- clm(CELI.Grade ~ cash_asset, data = daten_cash)
summary(modell_base_h4a)

#Full Model
modell_h4a <- clm(CELI.Grade ~ cash_asset + log_firm_size + industry_grouped + ROA + region,
                  data = daten_cash)
summary(modell_h4a)

run_diagnostics(modell_h4a, daten_cash, "H4a")
print(odds_ratios(modell_h4a))
ame_h4a <- avg_slopes(modell_h4a, variables = "cash_asset")
print(ame_h4a)

# H4b: Gearing (D/E) (tbd in Results Chapter 6)
daten_gear <- daten[!is.na(as.numeric(daten$Gearing.2021)), ]
daten_gear$gearing <- as.numeric(daten_gear$Gearing.2021)
cat("H4b Beobachtungen:", nrow(daten_gear), "von", nrow(daten), "\n")
cat("Coverage:", round(nrow(daten_gear)/nrow(daten)*100, 1), "%\n\n")
print(summary(daten_gear$gearing))

daten_gear <- prepare_controls(daten_gear)

#Base Model
modell_base_h4b <- clm(CELI.Grade ~ gearing, data = daten_gear)
summary(modell_base_h4b)

#Full Model
modell_h4b <- clm(CELI.Grade ~ gearing + log_firm_size + industry_grouped + ROA + region,
                  data = daten_gear)
summary(modell_h4b)

run_diagnostics(modell_h4b, daten_gear, "H4b")
print(odds_ratios(modell_h4b))
ame_h4b <- avg_slopes(modell_h4b, variables = "gearing")
print(ame_h4b)


# H4c: Current Ratio (tbd in Results Chapter 6)
daten_cr <- daten[!is.na(as.numeric(daten$Current.Ratio.2021)), ]
daten_cr$current_ratio <- as.numeric(daten_cr$Current.Ratio.2021)
cat("H4c Beobachtungen:", nrow(daten_cr), "von", nrow(daten), "\n")
cat("Coverage:", round(nrow(daten_cr)/nrow(daten)*100, 1), "%\n\n")
print(summary(daten_cr$current_ratio))

daten_cr <- prepare_controls(daten_cr)

#Base Model
modell_base_h4c <- clm(CELI.Grade ~ current_ratio, data = daten_cr)
summary(modell_base_h4c)

#Full Model
modell_h4c <- clm(CELI.Grade ~ current_ratio + log_firm_size + industry_grouped + ROA + region,
                  data = daten_cr)
summary(modell_h4c)

run_diagnostics(modell_h4c, daten_cr, "H4c")
print(odds_ratios(modell_h4c))
ame_h4c <- avg_slopes(modell_h4c, variables = "current_ratio")
print(ame_h4c)

#Zusammenfassung Tabelle

modelle <- list(H1a=modell_h1a, H1b=modell_h1b, H1c=modell_h1c, H2=modell_h2,
                H3a=modell_h3a, H3b=modell_h3b, H3c=modell_h3c,
                H4a=modell_h4a, H4b=modell_h4b, H4c=modell_h4c)

tab <- do.call(rbind, lapply(names(modelle), function(nm) {
  m   <- modelle[[nm]]
  b   <- names(m$beta)                       # Fokus-Prädiktoren
  est <- m$beta
  se  <- sqrt(diag(vcov(m)))[b]
  z   <- est / se
  p   <- 2 * pnorm(-abs(z))
  data.frame(
    Modell   = nm,
    Variable = b,
    OR       = round(exp(est), 3),
    p        = round(p, 3),
    N        = m$n,
    AIC      = round(AIC(m), 1),
    BIC      = round(BIC(m), 1),
    row.names = NULL
  )
}))
tab
fokus <- c("rev_share_pct","asset_share_pct","emp_share_pct",
           "sanction1","sanction2","log_subs","log_emp_ru","log_assets_ru",
           "cash_asset","gearing","current_ratio")
tab_fokus <- tab[tab$Variable %in% fokus, ]
tab_fokus

# ROBUSTNESS CHECK 1: 3-Level Collapsed Model (adressed in Section 5.3.2 )

daten$CELI.Grade_3 <- ifelse(daten$CELI.Grade %in% c("Digging In", "Buying Time"), "Stay",
                             ifelse(daten$CELI.Grade == "Scaling Back", "Partial", "Exit"))
daten$CELI.Grade_3 <- factor(daten$CELI.Grade_3,
                             levels = c("Stay", "Partial", "Exit"),
                             ordered = TRUE)

cat("Distribution 5-Level:\n")
print(table(daten$CELI.Grade))
cat("\nDistribution 3-Level:\n")
print(table(daten$CELI.Grade_3))
print(round(prop.table(table(daten$CELI.Grade_3)) * 100, 1))

# Subset: require valid subsidiaries AND sanction rating
daten_col <- daten[!is.na(as.numeric(daten$Subsidiaries_RU_Latest)) & 
                     !is.na(daten$Sanction_Rating), ]
daten_col$log_subs <- log(as.numeric(daten_col$Subsidiaries_RU_Latest) + 1)
daten_col$sanction <- factor(daten_col$Sanction_Rating)
daten_col$log_firm_size <- log(as.numeric(daten_col$Global_Revenue_2021_m) + 1)
daten_col$ROA <- as.numeric(daten_col$ROA.using.Profit...Loss..before.tax.2021....)
daten_col$region <- ifelse(daten_col$Country.CELI == "United States", "US",
                           ifelse(daten_col$Country.CELI == "United Kingdom", "UK",
                                  ifelse(daten_col$Country.CELI == "Germany", "Germany",
                                         ifelse(daten_col$Country.CELI == "France", "France",
                                                ifelse(daten_col$Country.CELI == "Finland", "Finland", "Other EU")))))
daten_col$region <- factor(daten_col$region)
daten_col$region <- relevel(daten_col$region, ref = "US")

daten_col <- daten_col[complete.cases(daten_col[, c("CELI.Grade_3", "log_subs", "sanction", 
                                                    "log_firm_size", "ROA", "region")]), ]
cat("Observations (complete cases):", nrow(daten_col), "\n\n")

#Model A: 3-Level with full region (may have separation for Finland)
cat("--- Model A: 3-Level with full region ---\n")
modell_3L <- clm(CELI.Grade_3 ~ log_subs + sanction + log_firm_size + ROA + region, data = daten_col)
summary(modell_3L)
print(nominal_test(modell_3L))

# Model B: 3-Level with binary region (US vs. Europe)
cat("\n--- Model B: 3-Level with binary region (primary robustness model) ---\n")
daten_col$region_bin <- ifelse(daten_col$Country.CELI == "United States", "US", "Europe")
daten_col$region_bin <- factor(daten_col$region_bin)
modell_3L_bin <- clm(CELI.Grade_3 ~ log_subs + sanction + log_firm_size + ROA + region_bin, data = daten_col)
summary(modell_3L_bin)

# Diagnostics for binary region model
print(nominal_test(modell_3L_bin))
cat("\n3-Level Distribution:\n")
print(table(daten_col$CELI.Grade_3))
print(round(prop.table(table(daten_col$CELI.Grade_3)) * 100, 1))

modell_null_3L <- clm(CELI.Grade_3 ~ 1, data = modell_3L_bin$model)
print(anova(modell_null_3L, modell_3L_bin))
ll_full <- as.numeric(logLik(modell_3L_bin))
ll_null <- as.numeric(logLik(modell_null_3L))
n <- nobs(modell_3L_bin)
nagelkerke_r2 <- (1 - exp((2/n) * (ll_null - ll_full))) / (1 - exp((2/n) * ll_null))
cat("Nagelkerke Pseudo R²:", round(nagelkerke_r2, 4), "\n")

print(odds_ratios(modell_3L_bin))
ame_3L_subs <- avg_slopes(modell_3L_bin, variables = "log_subs")
print(ame_3L_subs)
ame_3L_sanc <- avg_slopes(modell_3L_bin, variables = "sanction")
print(ame_3L_sanc)

# --- Coefficient comparison: 5-Level vs 3-Level (same data, same predictors) ---
cat("\n--- Coefficient comparison: 5-Level vs 3-Level ---\n")
modell_5L <- clm(CELI.Grade ~ log_subs + sanction + log_firm_size + ROA + region_bin, data = daten_col)
cat("5-Level coefficients:\n")
print(round(coef(modell_5L), 4))
cat("\n3-Level coefficients:\n")
print(round(coef(modell_3L_bin), 4))


# ROBUSTNESS CHECK 2: Combined Model (all significant predictors) As discussed in 5.3.2

#Subsidiaries + Sanctions (the two strongest predictors)
daten_comb <- daten[!is.na(as.numeric(daten$Subsidiaries_RU_Latest)) & 
                      !is.na(daten$Sanction_Rating), ]
daten_comb$log_subs <- log(as.numeric(daten_comb$Subsidiaries_RU_Latest) + 1)
daten_comb$sanction <- factor(daten_comb$Sanction_Rating)
daten_comb$log_firm_size <- log(as.numeric(daten_comb$Global_Revenue_2021_m) + 1)
daten_comb$ROA <- as.numeric(daten_comb$ROA.using.Profit...Loss..before.tax.2021....)
daten_comb$region <- ifelse(daten_comb$Country.CELI == "United States", "US",
                            ifelse(daten_comb$Country.CELI == "United Kingdom", "UK",
                                   ifelse(daten_comb$Country.CELI == "Germany", "Germany",
                                          ifelse(daten_comb$Country.CELI == "France", "France",
                                                 ifelse(daten_comb$Country.CELI == "Finland", "Finland", "Other EU")))))
daten_comb$region <- factor(daten_comb$region)
daten_comb$region <- relevel(daten_comb$region, ref = "US")

cat("--- Combined Model 1: log_subs + sanction ---\n")
modell_comb1 <- clm(CELI.Grade ~ log_subs + sanction + log_firm_size + ROA + region, data = daten_comb)
summary(modell_comb1)

#Diagnostics
print(nominal_test(modell_comb1))
modell_null_comb <- clm(CELI.Grade ~ 1, data = modell_comb1$model)
print(anova(modell_null_comb, modell_comb1))
ll_full <- as.numeric(logLik(modell_comb1))
ll_null <- as.numeric(logLik(modell_null_comb))
n <- nobs(modell_comb1)
nagelkerke_r2 <- (1 - exp((2/n) * (ll_null - ll_full))) / (1 - exp((2/n) * ll_null))
cat("Nagelkerke Pseudo R²:", round(nagelkerke_r2, 4), "\n")
print(odds_ratios(modell_comb1))
ame_comb1 <- avg_slopes(modell_comb1, variables = c("log_subs", "sanction"))
print(ame_comb1)

#AIC comparison: combined vs. individual models (same data)
modell_subs_only <- clm(CELI.Grade ~ log_subs + log_firm_size + ROA + region, data = modell_comb1$model)
modell_sanc_only <- clm(CELI.Grade ~ sanction + log_firm_size + ROA + region, data = modell_comb1$model)
cat("\nAIC Comparison:\n")
cat("Combined (subs + sanc):", AIC(modell_comb1), "\n")
cat("Subs only:            ", AIC(modell_subs_only), "\n")
cat("Sanction only:        ", AIC(modell_sanc_only), "\n")
cat("Null model:           ", AIC(modell_null_comb), "\n")

#add current_ratio
daten_comb$current_ratio <- as.numeric(daten_comb$Current.Ratio.2021)
daten_comb2 <- daten_comb[!is.na(daten_comb$current_ratio), ]

cat("\n--- Combined Model 2: + current_ratio ---\n")
modell_comb2 <- clm(CELI.Grade ~ log_subs + sanction + current_ratio + log_firm_size + ROA + region, 
                    data = daten_comb2)
summary(modell_comb2)

# LR test: Model 2 vs Model 1 (same data)
modell_comb1_refit <- clm(CELI.Grade ~ log_subs + sanction + log_firm_size + ROA + region, 
                          data = modell_comb2$model)
print(anova(modell_comb1_refit, modell_comb2))
cat("AIC Model 1 (refit):", AIC(modell_comb1_refit), "\n")
cat("AIC Model 2:        ", AIC(modell_comb2), "\n")

#add log_emp_ru + log_assets_ru (full model, smallest sample) ---
daten_comb$log_emp_ru <- log(as.numeric(daten_comb$Employees_Subsidiaries_RU) + 1)
daten_comb$log_assets_ru <- log(as.numeric(daten_comb$Assets_RU_Subsidiaries) + 1)
daten_comb3 <- daten_comb[!is.na(daten_comb$current_ratio) & 
                            !is.na(daten_comb$log_emp_ru) & 
                            !is.na(daten_comb$log_assets_ru), ]
cat("\n--- Combined Model 3: + log_emp_ru + log_assets_ru (full combined) ---\n")
cat("Observations:", nrow(daten_comb3[complete.cases(daten_comb3[, c("CELI.Grade", "log_subs", "sanction", 
                                                                     "current_ratio", "log_emp_ru", "log_assets_ru", "log_firm_size", "ROA", "region")]), ]), "\n")

modell_comb3 <- clm(CELI.Grade ~ log_subs + sanction + current_ratio + log_emp_ru + log_assets_ru + 
                      log_firm_size + ROA + region, data = daten_comb3)
summary(modell_comb3)

# VIF (measures how much each predictor's variance is inflated due to correlation with other predictors in the model)
library(car)
lm_helper <- lm(as.numeric(CELI.Grade) ~ log_subs + sanction + current_ratio + 
                  log_emp_ru + log_assets_ru + log_firm_size + ROA + region, data = daten_comb3)
cat("\nVIF (Variance Inflation Factors):\n")
print(vif(lm_helper))

# LR test: Model 3 vs Model 2 (same data)
modell_comb2_refit <- clm(CELI.Grade ~ log_subs + sanction + current_ratio + log_firm_size + ROA + region, 
                          data = modell_comb3$model)
print(anova(modell_comb2_refit, modell_comb3))

# AIC comparison all three
cat("\nAIC Comparison (stepwise):\n")
modell_comb1_refit3 <- clm(CELI.Grade ~ log_subs + sanction + log_firm_size + ROA + region, 
                           data = modell_comb3$model)
cat("Model 1 (subs + sanc):              ", AIC(modell_comb1_refit3), "\n")
cat("Model 2 (+ current_ratio):          ", AIC(modell_comb2_refit), "\n")
cat("Model 3 (+ emp_ru + assets_ru):     ", AIC(modell_comb3), "\n")

# Diagnostics full combined model
print(nominal_test(modell_comb3))
modell_null_comb3 <- clm(CELI.Grade ~ 1, data = modell_comb3$model)
print(anova(modell_null_comb3, modell_comb3))
ll_full <- as.numeric(logLik(modell_comb3))
ll_null <- as.numeric(logLik(modell_null_comb3))
n <- nobs(modell_comb3)
nagelkerke_r2 <- (1 - exp((2/n) * (ll_null - ll_full))) / (1 - exp((2/n) * ll_null))
cat("Nagelkerke Pseudo R²:", round(nagelkerke_r2, 4), "\n")
print(odds_ratios(modell_comb3))
ame_comb3 <- avg_slopes(modell_comb3, variables = c("log_subs", "sanction", "current_ratio", "log_emp_ru", "log_assets_ru"))
print(ame_comb3)


