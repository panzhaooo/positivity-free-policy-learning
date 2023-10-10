## Diabetes 130-US hospitals for years 1999-2008
## https://archive.ics.uci.edu/ml/datasets/Diabetes+130-US+hospitals+for+years+1999-2008

## We follow the pre-processing steps in Fairlearn
## https://github.com/fairlearn/talks/blob/main/2021_scipy_tutorial/preprocess.py


library(tidyverse)
library(mice)

data <- read_csv("diabetic_data.csv")

DATA <- data %>%
  rename(primary_diagnosis = diag_1) %>%
  mutate(
    # Create Outcome variables
    readmit = as.numeric(readmitted != "NO"),
    readmit_30days = as.numeric(readmitted == "<30"),
    # Replace missing values and re-code categories
    # Clean various medical codes
    age = case_match(
      age,
      c("[0-10)", "[10-20)", "[20-30)") ~ "30 years or younger",
      c("[30-40)", "[40-50)", "[50-60)") ~ "30-60 years",
      c("[60-70)", "[70-80)", "[80-90)", "[90-100)") ~ "Over 60 years"
    ),
    payer_code = case_match(payer_code, "?" ~ "Unknown", .default = payer_code),
    race = case_match(race, "?" ~ "Unknown", .default = race),
    admission_source_id = case_match(
      admission_source_id,
      c(1, 2, 3) ~ "Referral",
      7 ~ "Emergency",
      .default = "Other"
    ),
    discharge_disposition_id = case_match(
      discharge_disposition_id,
      1 ~ "Discharged to Home",
      .default = "Other"
    ),
    # Re-code Medical Specialties and Primary Diagnosis
    medical_specialty = case_match(
      medical_specialty,
      "?" ~ "Missing",
      "InternalMedicine" ~ "InternalMedicine",
      "Emergency/Trauma" ~ "Emergency/Trauma",
      "Family/GeneralPractice" ~ "Family/GeneralPractice",
      "Cardiology" ~ "Cardiology",
      "Surgery" ~ "Surgery",
      .default = "Other"
    ),
    primary_diagnosis = case_when(
      grepl("[7][1-3][0-9]", primary_diagnosis) ~ "Musculoskeltal Issues",
      grepl("250.*", primary_diagnosis) ~ "Diabetes",
      grepl("[4][6-9][0-9]|[5][0-1][0-9]|786", primary_diagnosis) ~ "Respitory Issues",
      grepl("[5][8-9][0-9]|[6][0-2][0-9]|788", primary_diagnosis) ~ "Genitourinary Issues",
      .default = "Other"
    ),
    #Binarize and bin features
    medicare = as.numeric(payer_code == "MC"),
    medicaid = as.numeric(payer_code == "MD"),
    had_emergency = as.numeric(number_emergency > 0),
    had_inpatient_days = as.numeric(number_inpatient > 0),
    had_outpatient_days = as.numeric(number_outpatient > 0)) %>%
  select(race, gender, age, discharge_disposition_id, admission_source_id, 
         time_in_hospital, medical_specialty, num_lab_procedures, num_procedures,
         num_medications, primary_diagnosis, number_diagnoses, max_glu_serum,
         A1Cresult, insulin, change, diabetesMed, medicare, medicaid, 
         had_emergency, had_inpatient_days, had_outpatient_days, readmitted,
         readmit, readmit_30days)


DATA <- DATA %>%
  mutate(
    race = case_match(
      race,
      "Caucasian" ~ 1,
      c("AfricanAmerican", "Asian", "Hispanic", "Other") ~ 0,
      "Unknown" ~ NA),
    gender = case_match(
      gender,
      "Female" ~ 1,
      "Male" ~ 0,
      "Unknown/Invalid" ~ NA
    ),
    age = case_match(
      age,
      c("30 years or younger", "30-60 years") ~ 1,
      "Over 60 years" ~ 0
    ),
    discharge_home = case_match(
      discharge_disposition_id,
      "Discharged to Home" ~ 1,
      "Other" ~ 0
    ),
    admission_emergency = case_match(
      admission_source_id,
      "Emergency" ~ 1,
      .default = 0
    ),
    admission_referral = case_match(
      admission_source_id,
      "Referral" ~ 1,
      .default = 0
    ),
    insulin_steady = case_match(
      insulin,
      "Steady" ~ 1,
      .default = 0
    ),
    change = case_match(
      change,
      "Ch" ~ 1,
      "No" ~ 0
    ),
    diabetesMed = case_match(
      diabetesMed,
      "Yes" ~ 1,
      "No" ~ 0
    )) %>%
  select(race, gender, age, discharge_home, admission_emergency, admission_referral,
         time_in_hospital, num_lab_procedures, num_procedures,
         num_medications, number_diagnoses,
         insulin_steady, change, diabetesMed, medicare, medicaid, 
         had_emergency, had_inpatient_days, had_outpatient_days, readmitted,
         readmit, readmit_30days)


set.seed(9907)
imp <- mice(DATA)
DATA <- complete(imp)

save(DATA, file = "Diabetes.RData")


