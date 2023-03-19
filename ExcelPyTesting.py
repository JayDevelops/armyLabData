import openpyxl
import pandas as pd

book = openpyxl.load_workbook(r'Copy of ADM_-_from_Joel_-_Sept-2013.xlsx')
sheet = book["Model"]

print(sheet["A8"].value)

for row in sheet.iter_rows(min_row=8, max_row=31, min_col=1, max_col=1, values_only=True):
    print(row)

#adm_model_df = pd.read_excel('Copy of ADM_-_from_Joel_-_Sept-2013.xlsx', sheet_name='Model')
#adm_model_freqHz = pd.read_excel('Copy of ADM_-_from_Joel_-_Sept-2013.xlsx', skiprows = 7, nrows= 23, usecols = 'A', sheet_name='Model')
#print(adm_model_freqHz)