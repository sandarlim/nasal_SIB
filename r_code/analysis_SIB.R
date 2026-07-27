# Protocol: SIB vs regular
# =============

# Basic log-rank tests
# --------
survdiff(Surv(PFS, Progression_or_Death_y_n) ~ Protocol, data = dat)
survdiff(Surv(OS, Dead_y_n) ~ Protocol, data = dat)

# KM curves
# --------
KM_P_PFS <- survfit(Surv(PFS, Progression_or_Death_y_n) ~ Protocol, data = dat)
plot_KM_P_PFS <- ggsurvplot(KM_P_PFS,
                             pval = TRUE, conf.int = TRUE,
                             risk.table = TRUE, # Add risk table
                             risk.table.col = "strata", # Change risk table color by groups
                             linetype = "strata", # Change line type by groups
                             surv.median.line = "hv", # Specify median survival
                             ggtheme = theme_bw(), # Change ggplot2 theme
                             palette = c("#868686FF", "#0073C2FF"),
                             legend.labs = c("regular protocol", "boost protocol"),
                             ylab = "Survival probability (PFS)")
plot_KM_P_PFS # 700-500

KM_P_OS <- survfit(Surv(OS, Dead_y_n) ~ Protocol, data = dat)
ggsurvplot(KM_P_OS,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#868686FF", "#0073C2FF"),
           legend.labs = c("regular protocol", "boost protocol"),
           ylab = "Survival probability (OS)")


# Stage
# =============
# Basic log-rank tests
# --------
survdiff(Surv(PFS, Progression_or_Death_y_n) ~ Stage4_y_n, data = dat)
survdiff(Surv(OS, Dead_y_n) ~ Stage4_y_n, data = dat)

KM_stage_PFS <- survfit(Surv(PFS, Progression_or_Death_y_n) ~ Stage4_y_n, data = dat)
plot_KM_stage_PFS <- ggsurvplot(KM_stage_PFS,
                                 pval = TRUE, conf.int = TRUE,
                                 risk.table = TRUE, # Add risk table
                                 risk.table.col = "strata", # Change risk table color by groups
                                 linetype = "strata", # Change line type by groups
                                 surv.median.line = "hv", # Specify median survival
                                 ggtheme = theme_bw(), # Change ggplot2 theme
                                 palette = c("#868686FF", "#0073C2FF"),
                                 legend.labs = c("Stage I-III", "Stage IV"),
                                 ylab = "Survival probability (PFS)")
plot_KM_stage_PFS

KM_stage_OS <- survfit(Surv(OS, Dead_y_n) ~ Stage4_y_n, data = dat)
plot_KM_stage_OS <- ggsurvplot(KM_stage_OS,
                                pval = TRUE, conf.int = TRUE,
                                risk.table = TRUE, # Add risk table
                                risk.table.col = "strata", # Change risk table color by groups
                                linetype = "strata", # Change line type by groups
                                surv.median.line = "hv", # Specify median survival
                                ggtheme = theme_bw(), # Change ggplot2 theme
                                palette = c("#868686FF", "#0073C2FF"),
                                legend.labs = c("Stage I-III", "Stage IV"),
                                ylab = "Survival probability (OS)")
plot_KM_stage_OS


# Univariate analysis of progression
# ======================
# PFS
p <- coxph(Surv(PFS, Progression_or_Death_y_n) ~ Protocol, data = dat)
summary(p)
age <- coxph(Surv(PFS, Progression_or_Death_y_n) ~ Age + strata(Protocol), data = dat)
summary(age)
weight <- coxph(Surv(PFS, Progression_or_Death_y_n) ~ Weight + strata(Protocol), data = dat)
summary(weight)
stage4 <- coxph(Surv(PFS, Progression_or_Death_y_n) ~ Stage4_y_n + strata(Protocol), data = dat)
summary(stage4)
epistaxis <- coxph(Surv(PFS, Progression_or_Death_y_n) ~ Epistaxis + strata(Protocol), data = dat)
summary(epistaxis)
tumor_diff <- coxph(Surv(PFS, Progression_or_Death_y_n) ~ Tumor_undifferentiated + strata(Protocol), data = dat)
summary(tumor_diff)
tumor_code <- coxph(Surv(PFS, Progression_or_Death_y_n) ~ TumorCode + strata(Protocol), data = dat)
summary(tumor_code)
tumor_vol <- coxph(Surv(PFS, Progression_or_Death_y_n) ~ GTV + strata(Protocol), data = dat)
summary(tumor_vol)

# OS
p2 <- coxph(Surv(OS, Dead_y_n) ~ Protocol, data = dat)
summary(p2)
age2 <- coxph(Surv(OS, Dead_y_n) ~ Age + strata(Protocol), data = dat)
summary(age2)
weight2 <- coxph(Surv(OS, Dead_y_n) ~ Weight + strata(Protocol), data = dat)
summary(weight2)
stage42 <- coxph(Surv(OS, Dead_y_n) ~ Stage4_y_n + strata(Protocol), data = dat)
summary(stage42)
epistaxis2 <- coxph(Surv(OS, Dead_y_n) ~ Epistaxis + strata(Protocol), data = dat)
summary(epistaxis2)
tumor_diff2 <- coxph(Surv(OS, Dead_y_n) ~ Tumor_undifferentiated + strata(Protocol), data = dat)
summary(tumor_diff2)
tumor_code2 <- coxph(Surv(OS, Dead_y_n) ~ TumorCode + strata(Protocol), data = dat)
summary(tumor_code2)
tumor_vol2 <- coxph(Surv(OS, Dead_y_n) ~ GTV + strata(Protocol), data = dat)
summary(tumor_vol2)


# Full Cox model
# ------------------
table(dat$Progression_or_Death_y_n) # 34 events
m_pfs0 <- coxph(Surv(PFS, Progression_or_Death_y_n) ~
                  Protocol +
                  Stage4_y_n +
                  Age +
                  Weight +
                  Epistaxis +
                  GTV,
                data = dat)
summary(m_pfs0)
summary(m_pfs0)$coefficients

table(dat$Dead_y_n) # 39 events
m_os0 <- coxph(Surv(OS, Dead_y_n) ~ Protocol +
                 Stage4_y_n +
                 Age +
                 Weight +
                 Epistaxis +
                 GTV,
               data = dat)
summary(m_os0)
summary(m_os0)$coefficients
