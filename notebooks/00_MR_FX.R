#### SMR & HEIDI Functions

### SMR p-value given the effects of variants on exposure x and outcome y
pval_smr <- function(b_zx, se_zx, b_zy, se_zy){
    ## Z - scores squared
    z2_zy <- (b_zy / se_zy)^2
    z2_zx <- (b_zx / se_zx)^2
    ## T statistic
    t <- (z2_zy * z2_zx) / (z2_zy + z2_zx)
    ## P value based on Chi squared distribution
    ## 1 degree of freedom
    p <- pchisq(t, df = 1, lower.tail = F)
    return(p)
}

### SMR standard error of beta
std_err_smr <- function(b_zx, se_zx, b_zy, se_zy){
    ## SMR beta
    b_xy <- b_zy / b_zx
    ## Squared coefficients of variation
    varcoef_zx <- (se_zx/b_zx)^2
    varcoef_zy <- (se_zy/b_zy)^2
    ## Sum of coefficients of variation
    sum_varcoef <- varcoef_zx + varcoef_zy
    ## Variance
    var_xy <- (b_xy^2) * sum_varcoef
    se <- sqrt(var_xy)
    return(se)
}

### Covariance of bXY
cov_bXY <- function(bzx, bzx_se, bzy, bzy_se, ldrho) {
    bXY = bzy / bzx
    zscoreZX = bzx / bzx_se
    zszxij = zscoreZX %*% t(zscoreZX)
    bxyij = bXY %*% t(bXY)
    sezyij = bzy_se %*% t(bzy_se)
    bzxij = bzx %*% t(bzx)
    (ldrho*sezyij/bzxij) + (ldrho*bxyij/zszxij) - (bxyij/(zszxij^2))
}

### HEIDI p-value
heidi_pvalue <- function(bzx, bzx_se, bzy, bzy_se, ldrho, topsnp_index) {
    ## 1. Top SNP
    bxy_top <- bzy[topsnp_index]/bzx[topsnp_index]
    bxy_top_se <- std_err_smr(bzx[topsnp_index], bzx_se[topsnp_index],
                              bzy[topsnp_index], bzy_se[topsnp_index])
    ## 2. d = bxy_i - bxy_top
    bxy_snp <- bzy / bzx
    d <- bxy_snp - bxy_top
    cov_bxy <- cov_bXY(bzx, bzx_se, bzy, bzy_se, ldrho)
    var_d <- cov_bxy - cov_bxy[topsnp_index,] + bxy_top_se^2
    var_d <- t(t(var_d) - cov_bxy[topsnp_index,])
    ## 3. Removing top SNP
    d <- d[-topsnp_index]
    var_d <- var_d[-topsnp_index, -topsnp_index]
    m <- length(bzx) - 1
    ## 4. Chi-square values - Zd
    chival <- d^2/diag(var_d)
    ## 5. R - Correlation matrix of Zd
    corr_d <- diag(m)
    for( i in 1 : (m-1) ) {
        for( j in (i+1) : m ) {
            corr_d[i,j] = corr_d[j,i] = var_d[i,j] / sqrt(var_d[i,i]*var_d[j,j])
        }         
    }
    ## 6. Eigen decomposition
    lambda <- eigen(corr_d, symmetric=TRUE, only.values=TRUE)$values
    ## 7. Estimate p-value using saddlepoint method
    t <- length(lambda)
    survey::pchisqsum(sum(chival), df=rep(1,t), a=lambda,
                      method="saddlepoint", lower.tail = F)
}