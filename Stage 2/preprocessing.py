
import pandas as pd
import numpy as np

# Load the dataset
file_path = 'path_to_your_file.csv'  # Update with your file path
data = pd.read_csv(file_path, sep='\t')

# Step 1: Whitespace and Character Cleanup
# Strip leading/trailing whitespace from column names
data.columns = data.columns.str.strip()

# Strip leading/trailing whitespace from data entries
data = data.applymap(lambda x: x.strip() if isinstance(x, str) else x)

# Step 2: Convert Different Forms of NAs to Standard 'NA'
# Replace common NA values with np.nan for consistency
na_values = ['-', '?', 'â€“', 'N/A', 'n/a', 'NA']
data.replace(na_values, np.nan, inplace=True)

# Step 3: Remove Irrelevant or Redundant Columns
# Based on analysis, remove columns that are not useful
columns_to_remove = ['Alternative name', 'Clinical trials']  # Add columns as identified
data.drop(columns=columns_to_remove, axis=1, inplace=True)

# Step 4: Handle Missing Data
# Option 1: Replace '-' and '?' with np.nan (already done in Step 2)
# Option 2: Remove rows with excessive missing values (threshold-based removal)

# Set a threshold for missing values: if more than 50% of a row is missing, drop it
threshold = 0.5 * data.shape[1]  # 50% of total columns
data.dropna(thresh=threshold, axis=0, inplace=True)

# Step 5: Further Clean Specific Characters
# Replace instances of â€“ or other special characters with np.nan
data.replace({'â€“': np.nan}, inplace=True)

# Step 6: Fix Formatting of Key Columns (if needed)
# Example: Ensure numeric columns are formatted correctly or convert types if necessary
# Example: data['SomeNumericColumn'] = pd.to_numeric(data['SomeNumericColumn'], errors='coerce')

# Step 7: Final Review of Data Types (optional)
# Review the data types of all columns and ensure consistency
data.info()

# Save the cleaned dataset
cleaned_file_path = 'path_to_save_cleaned_data.csv'
data.to_csv(cleaned_file_path, index=False)

print(f"Data preprocessing completed and saved to {cleaned_file_path}")
