### Plotting

#### power bar
power_bar <- function(effect_size, res_list, cutoff=0.05,
                      eff_gap=c(0.5, 1.5, 2.5, 3.5, 4.5, 6.5),
                      Methods = c("dcats_bin_Wald", "dcats_BB_Wald",
                                  "std_BB_LRT", "dcats_BB_LRT")) {
    group_idx = c()
    for (eff_tmp in effect_size) {
        group_idx <- c(group_idx, which.min(abs(eff_gap - abs(eff_tmp))))
    }

    TPR_list <- method_list <- effsize_list <- c()
    for (i in seq_len(length(eff_gap))) {
        idx <- which(group_idx == i)
        print(paste(i, idx))
        for (method_tmp in Methods) {
            TPR_list <- c(TPR_list, mean(res_list[[method_tmp]][, idx] < cutoff))
            method_list <- c(method_list, method_tmp)
            effsize_list <- c(effsize_list, eff_gap[i])
        }
    }
    df.fig = data.frame(effect_size_abs = factor(effsize_list, levels=eff_gap),
                        method = factor(method_list, levels = Methods),
                        power=TPR_list)
    ggplot(df.fig, aes(fill=method, y=power, x=effect_size_abs)) +
        geom_bar(position="dodge", stat="identity") +
        viridis::scale_fill_viridis(discrete=TRUE) + theme_bw()
}

##### qq plot
qq_plot <- function(pvalues) {
    obs_pval <- as.numeric(pvalues)
    exp_pval <- (rank(obs_pval, ties.method="first")+.5)/(length(obs_pval)+1)
    df.fig <- data.frame(obs_pval = obs_pval, exp_pval = exp_pval)
    ggplot(data = df.fig, aes(x=-log10(exp_pval), y=-log10(obs_pval))) +
        geom_point(shape=1) +
        geom_abline(intercept = 0, slope = 1, color="grey") +
        theme_bw()
}

FPR_qqplot <- function(res_list,
                       method_list = c("dcats_bin_LRT",
                                       "dcats_BB_Wald",
                                       "std_BB_LRT",
                                       "dcats_BB_LRT")) {
    pp_list <- list()
    for (i in seq_len(length(method_list))) {
        pp_list[[i]] <- qq_plot(res_list[[method_list[i]]]) +
            ggtitle(method_list[i])
    }
    ggpubr::ggarrange(pp_list[[1]], pp_list[[2]], pp_list[[3]], pp_list[[4]],
                      nrow=1, ncol = 4)
}





### Main simulator
sys_simulator <- function(prop1, prop2, simMM, n_sample=50,
                          concertration1=50, concertration2=50,
                          n_rep1=4, n_rep2=4,
                          n_cell1=2500, n_cell2=2500, cv=0.2,
                          keep_betabin=TRUE) {
    mat_dcats_old <- mat_dcats_bin <- mat_bin_LRT <-
        mat_std_bb <- mat_bb_wald <- mat_bb_LRT <-
        matrix(NA, n_sample, length(prop1))

    for (ir in seq_len(n_sample)) {
        print(paste("iteration", ir))
        total_samples1 = round(rnorm(n_rep1, n_cell1, n_cell1 * cv))
        total_samples2 = round(rnorm(n_rep2, n_cell2, n_cell2 * cv))

        dirichlet_s1 = prop1 * concertration1
        dirichlet_s2 = prop2 * concertration2

        sim_dat <- DCATS::simulator_base(total_samples1, total_samples2,
                                         dirichlet_s1, dirichlet_s2, simMM)

        res_old <- dcats_fit(sim_dat$numb_cond1, sim_dat$numb_cond2, simMM)
        res_bin <- dcats_betabin(sim_dat$numb_cond1, sim_dat$numb_cond2,
                                 simMM, binom_only = TRUE)
        mat_dcats_old[ir, ] <- res_old$pvals
        mat_dcats_bin[ir, ] <- res_bin$pvals
        mat_bin_LRT[ir, ] <- res_bin$LRT_pvals

        if (keep_betabin) {
            res_std_bb <- betabinLRT(sim_dat$numb_cond1, sim_dat$numb_cond2)
            mat_std_bb[ir, ] <- res_std_bb$pvals

            res_dcats_bb <- dcats_betabin(sim_dat$numb_cond1, sim_dat$numb_cond2,
                                        simMM, binom_only = FALSE)
            mat_bb_wald[ir, ] <- res_dcats_bb$pvals
            mat_bb_LRT[ir, ]  <- res_dcats_bb$LRT_pvals
        }
        # res.dcatbb_uni <- dcats_betabin(sim_dat$numb_cond1, sim_dat$numb_cond1,
        #                                 simMM_uni, binom_only = TRUE)
    }
    list("std_BB_LRT"=mat_std_bb,
         "bin_Wald_old"=mat_dcats_old,
         "dcats_bin_Wald"=mat_dcats_bin,
         "dcats_bin_LRT"=mat_bin_LRT,
         "dcats_BB_Wald"=mat_bb_wald,
         "dcats_BB_LRT"=mat_bb_LRT)
}


library(DCATS)
library(ggpubr)
count_mat <- as.matrix(read.table("examples/ALM_data/subjectinfomat_AML.txt"))
simMM <- as.matrix(read.table("examples/ALM_data/numemisclassM_AML.txt"))
simMM_uni <- get_similarity_mat(K=21, confuse_rate = 0.05)

idx_row <- idx_col <- match(colnames(count_mat), colnames(simMM))
simMM <- simMM[, idx_col][idx_row, ]


prop1 <- (count_mat[1, ] + 1) / sum(count_mat[1, ] + 1)
prop2 <- (count_mat[5, ] + 1) / sum(count_mat[5, ] + 1)

n_rep1 <- n_rep2 <- 10
n_cell1 <- n_cell2 <- 2500
concertration1 <- concertration2 <- 50



### False positive controls
res_FP_ctrl1 <- sys_simulator(prop2, prop2, simMM, n_sample = 10,
                              n_rep1 = 10, n_rep2=10)
res_FP_ctrl2 <- sys_simulator(prop2, prop2, simMM, n_sample = 10,
                              n_rep1 = 3, n_rep2=3)

FPR_qqplot(res_FP_ctrl1)
FPR_qqplot(res_FP_ctrl2)


### Power analysis
eff_size = qlogis(prop1) - qlogis(prop2) # logit <- qlogis
res_power1 <- sys_simulator(prop1, prop2, simMM, n_sample=10,
                            n_rep1 = 10, n_rep2=10)
res_power2 <- sys_simulator(prop1, prop2, simMM, n_sample=10,
                            n_rep1 = 3, n_rep2=3)

power_bar(eff_size, res_power1)
power_bar(eff_size, res_power2)




