library(tidyverse)
library(cbioportalR)

set.seed(20240101)

set_cbioportal_db("public")

# get clinical data for the brca_metabric study
patients <- available_patients("brca_metabric")
clinical <- get_clinical_by_patient(
    patient_id = patients$patientId, study_id = "brca_metabric"
)

# wrangle data
clinical <- clinical |> 
  select(-uniquePatientKey, -studyId) |> 
  pivot_wider(
    names_from = clinicalAttributeId, values_from = value,
  ) |> 
  mutate(across(everything(), parse_guess)) 

# get the mutation data. For convenience we'll limit the data to only genenames
# and protein change.
mutations <- get_mutations_by_study(study_id = "brca_metabric") |> 
  select(patientId, hugoGeneSymbol, proteinChange)

# prepare outcomes
outcomes <- clinical |> 
  select(patientId, OS_MONTHS, VITAL_STATUS) |> 
  arrange(patientId) |> 
  mutate(
    time  = (OS_MONTHS * 30.43) |> as.integer(),
    event = VITAL_STATUS |> 
      case_match(
        "Living" ~ 0L, "Died of Disease" ~ 1L, "Died of Other Causes" ~ 2L
      ),
    event_str = VITAL_STATUS, 
    .keep = "unused"
  ) |> 
  filter(!is.na(event), !is.na(time)) |> 
  mutate(
    split = sample(
      c("train", "test"), size = n(), prob = c(0.8, 1), replace = TRUE
    )
  )

## export tsv files

clinical |> 
  select(-SEX, -COHORT, -INTCLUST, -SAMPLE_COUNT, -SEX,
         -RFS_STATUS, -OS_STATUS, -VITAL_STATUS, -OS_MONTHS,
         -RFS_MONTHS) |> 
  semi_join(outcomes, join_by(patientId)) |> 
  write_tsv("metabric_clinical.tsv")

outcomes |> 
  write_tsv("metabric_survival.tsv")

mutations |> 
  semi_join(outcomes, join_by(patientId)) |> 
  write_tsv("metabric_mutations.tsv")