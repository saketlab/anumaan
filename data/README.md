# anumaan Data Directory

This directory contains the datasets required for AMR (Antimicrobial Resistance) burden analysis.

## Important Note

**Data files are not included in the Git repository** due to their size. You need to obtain and place the datasets in this folder before running the analysis.

## Required Datasets

Place the following datasets in this directory:

```
anumaan/data/
├── .gitkeep
├── README.md (this file)
└── [your data files here]
```

### Expected File Formats

The analysis pipeline supports the following data formats:
- CSV files (`.csv`)
- Excel files (`.xlsx`, `.xls`)
- R data files (`.rds`)
- Text files (`.txt`)

## How to Obtain the Data

<!-- TODO: Update this section with specific instructions for your datasets -->

### Option 1: Contact the Repository Maintainer

Contact the repository maintainer or principal investigator to obtain access to the datasets.

### Option 2: Download from External Source

If datasets are hosted externally (e.g., Zenodo, Google Drive, institutional repository):

1. Download the datasets from: [ADD LINK HERE]
2. Extract/place the files in this `anumaan/data/` directory
3. Verify the file names match those expected by the analysis scripts

### Option 3: Use Your Own Data

If you're using your own AMR surveillance data, ensure it follows the expected schema:

- **Required columns**: [TODO: List required columns]
- **Date format**: [TODO: Specify date format]
- **Coding standards**: [TODO: ICD-10, pathogen names, antibiotic codes, etc.]

## Data Structure

<!-- TODO: Document your specific data structure -->

Example expected structure:
```
- Patient ID
- Date of culture
- Pathogen identified
- Antibiotic tested
- Susceptibility result (S/I/R)
- Specimen type
- Demographics (age, gender, etc.)
```

## Verification

After placing your data files, verify the setup by running:

```r
# Load the package
library(anumaan)

# Check if data files are present
list.files("data/")

# Or from the package root
list.files("anumaan/data/")
```

## Data Privacy and Ethics

- **Do not commit** actual patient data to version control
- Ensure you have appropriate ethical approval and data sharing agreements
- Anonymize/de-identify data according to your institutional requirements
- Follow local data protection regulations (GDPR, HIPAA, etc.)

## Questions?

If you encounter issues with data access or format, please:
- Open an issue on the GitHub repository
- Contact: [ADD CONTACT EMAIL]
- Refer to the main README for additional documentation

---

*Last updated: January 2026*
