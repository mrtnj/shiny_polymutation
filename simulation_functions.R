
## Simulate variation
sim_variation <- function(N, mu, s, h, gen, loci, progress_function = NULL) {
  ## Matrix to store allele frequencies over generations
  q <- matrix(0, ncol = loci, nrow = gen)
  mean_fitness <- numeric(gen)
  mean_fitness[1] <- 1
  n_ind <- numeric(gen)
  n_ind[1] <- N
  for (i in 2:gen) {

    if (!is.null(progress_function)) {
        progress_function(1/gen)
    }

    geno <- matrix(0, ncol = loci, nrow = N)
    for (j in 1:loci) {
        if (q[i - 1, j] == 0) {
          geno[, j] <- 0
        } else {
          geno[, j] <- rbinom(N, size = 2, prob = q[i - 1, j])
        }
    }
    ## Calculate fitness as one factor for heterozygous loci and one for homozygous loci
    ## Keep only individuals that draw a random number above their fitness
    fitness <- (1 - h * s) ^ rowSums(geno == 1) * (1 - s) ^ rowSums(geno == 2)
    survival_roll <- runif(N, 0, 1)
    survival <- ifelse(survival_roll < fitness, TRUE, FALSE)
    surviving <- geno[survival,]

    ## Check if population has died
    if (nrow(surviving) < 2) {
        break
    }

    ## Mutate the genotypes of survivors
    mutation_roll <- matrix(runif(nrow(surviving) * loci, 0, 1),
                            nrow = nrow(surviving), ncol = loci)
    mutated <- surviving
    mutated[which(mutation_roll < mu & surviving == 0)] <- 1
    mutated[which(mutation_roll < mu & surviving == 1)] <- 2
    mutation_roll2 <- matrix(runif(nrow(surviving) * loci, 0, 1),
                             nrow = nrow(surviving), ncol = loci)
    mutated[which(mutation_roll2 < mu & mutated == 0)] <- 1
    mutated[which(mutation_roll2 < mu & mutated == 1)] <- 2

    ## Allele frequencies
    n_P = apply(mutated, 2, function(x) sum(x == 0))
    n_H = apply(mutated, 2, function(x) sum(x == 1))
    n_Q = apply(mutated, 2, function(x) sum(x == 2))
    q[i,] <- (n_H + n_Q * 2) / (2 * (n_P + n_H + n_Q))

    mean_fitness[i] <- mean(fitness)
    n_ind[i] <- nrow(surviving)
  }
  list(q = q, mean_fitness = mean_fitness, N = n_ind)
}

## Plot simulations
plot_simulations <- function(sim) {
  melted_q <- melt(sim$q)
  colnames(melted_q) <- c("gen", "locus", "q")
  ymax <- ifelse(max(melted_q$q) > 0.25, 1, 0.25)
  plot_q <- qplot(x = gen,
                  y = q,
                  group = locus,
                  data = melted_q,
                  geom = "line",
                  alpha = I(0.5)) +
    theme_bw() +
    ylim(0, ymax) +
    xlab("Generation") +
    ylab("Variant frequency")

  plot_fitness <- qplot(x = 1:length(sim$mean_fitness),
                        y = sim$mean_fitness,
                        geom = "line") +
    theme_bw() +
    xlab("Generation") +
    ylab("Average fitness")

  list(plot_fitness, plot_q)
}
