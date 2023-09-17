write.table(t(c(1, 2, 3)),
    file = "results/test.csv", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
)
cat("\n", file = "results/test.csv", append = TRUE)
write.table(t(c(1, 2, 3)),
    file = "results/test.csv", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
)
