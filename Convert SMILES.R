library(here)
library(rcdk)

SMILES_path <- here("Input","SMILES Code List Arzneistoffe und Metaboliten.txt")
Smiles_Arzneimittel <- read.csv(SMILES_path, header = FALSE)

smiles_vector <- Smiles_Arzneimittel $V1

convert_to_canonical <- function(smiles, index = NA) {
  mol <- tryCatch(parse.smiles(smiles)[[1]], error = function(e) NULL)
  if (!is.null(mol)) {
    get.smiles(mol, smiles.flavors("Canonical"))
  } else {
    cat("Fehlerhafte SMILES in Zeile", index, ":", smiles, "\n")
    return(NA)
  }
}

canonical_smiles <- mapply(convert_to_canonical, smiles_vector, seq_along(smiles_vector))

Smiles_Arzneimittel$canonical_smiles <- canonical_smiles

write_csv(
  data.frame(canonical_smiles = Smiles_Arzneimittel$canonical_smiles),
  here::here("Output", "canonical_SMILES_Arzneimittel_output.csv")
)
