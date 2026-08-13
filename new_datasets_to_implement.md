
original excel file is here: "/Users/christianhalter/Desktop/ECODA_PAPER_DATASETS.xlsx"

i guess this was the code the model used to read and extract the data?

```python
import openpyxl
import pandas as pd
import json

wb = openpyxl.load_workbook('ECODA_PAPER_DATASETS.xlsx', data_only=True)
print("Sheet names:", wb.sheetnames)

for sheetname in wb.sheetnames:
    ws = wb[sheetname]
    print(f"\n--- Sheet: {sheetname} ---")
    print("Dimensions:", ws.dimensions)
    
    # Read rows and check fill colors
    rows_data = []
    for row in ws.iter_rows(values_only=False):
        row_vals = []
        row_colors = []
        for cell in row:
            val = cell.value
            fill = cell.fill
            fg_color = None
            if fill and fill.fill_type:
                if fill.start_color:
                    fg_color = fill.start_color.rgb or fill.start_color.theme
            row_vals.append(val)
            row_colors.append(fg_color)
        if any(v is not None for v in row_vals):
            rows_data.append({'values': row_vals, 'colors': row_colors})
            
    print(f"Total non-empty rows: {len(rows_data)}")
    for i, r in enumerate(rows_data[:10]):
        print(f"Row {i}: {r['values'][:8]}")


```

```text
Sheet names: ['Sheet1']

--- Sheet: Sheet1 ---
Dimensions: A1:K13
Total non-empty rows: 13
Row 0: ['Benchmark', 'Batch effect', 'Data', 'No. cells', 'No. samples', 'No. cell_types', 'No. classes and Severity Scores', 'Possible exclusion criteria']
Row 1: ['X', None, 'Alzheimer [9]', '1 395 601', 83, 18, '2; (Control=1, Dementia=2)', None]
Row 2: [None, 'X (10x v2 and v3)\nSample prep', 'Breast cancer [19]', '714 331', 126, 10, '2; (Control=1, Cancer=2)', None]
Row 3: [None, 'X (collection site)', 'Covid-19 PBMC [5]', '993 171', 151, 10, '3; (Normal=1,Mild=2, Isch=3)', None]
Row 4: [None, 'X (different studies)\n10x v2 and v3', 'Diabetes [28]', '264 235', 52, 13, '4; (Normal=1, EPD =2, Types I and II =3)', 'Possible but:\n- 9 dataset\n- Lots of conditions (Diabetes, chemical stress, strain, etc.)']
Row 5: ['X (subset kidney, scRNA-seq, by disease)', 'X (need to check, according to reviewer multiple sites)', 'Kidney [7]', '104 314', 45, 14, '3; (Normal=1,CKD=2, AKI=3)', None]
Row 6: [None, None, 'Kidney cancer [29]', '50 236', 17, 6, '2; (Normal=1, Cancer=2)', 'Too few samples']
Row 7: ['X', None, 'Lupus PBMC [30]', '1 263 676', 261, 11, '2; (Normal=1, SLE=2)', None]
Row 8: [None, None, 'Lung [4]', '941 504', 165, 12, '5; (Normal=1, COPD=2, SCLC, LA, NSCLC=3)', '"- Too many technical conditions (as described in paper)"']
Row 9: ['X', None, 'Myocardial infraction (2) [6]', '132 888', 23, 11, '3; (Normal=1, Isch=2, Fib=3)', 'Many bio conditions (location + disease)']


```

```python
import openpyxl

wb = openpyxl.load_workbook('ECODA_PAPER_DATASETS.xlsx', data_only=True)
ws = wb['Sheet1']

data = []
for r_idx, row in enumerate(ws.iter_rows(values_only=False), start=1):
    row_info = []
    for c_idx, cell in enumerate(row, start=1):
        color_hex = None
        if cell.fill and cell.fill.fill_type:
            # openpyxl start_color can be RGB hex or ARGB
            color_hex = cell.fill.start_color.rgb
        row_info.append({
            'col': c_idx,
            'val': cell.value,
            'color': color_hex
        })
    data.append(row_info)

for r_idx, row in enumerate(data):
    vals = [f"C{c['col']}: {repr(c['val'])} (col:{c['color']})" for c in row if c['val'] is not None or c['color'] is not None]
    print(f"Row {r_idx}:")
    for v in vals:
        print("  ", v)


```

```text
Row 0:
   C1: 'Benchmark' (col:None)
   C2: 'Batch effect' (col:None)
   C3: 'Data' (col:None)
   C4: 'No. cells' (col:None)
   C5: 'No. samples' (col:None)
   C6: 'No. cell_types' (col:None)
   C7: 'No. classes and Severity Scores' (col:None)
   C8: 'Possible exclusion criteria' (col:None)
   C9: 'Comment' (col:None)
   C10: 'Has high-granularity cell types' (col:None)
   C11: 'Available on cellxgene' (col:None)
Row 1:
   C1: 'X' (col:None)
   C3: 'Alzheimer [9]' (col:Values must be of type <class 'str'>)
   C4: '1 395 601' (col:Values must be of type <class 'str'>)
   C5: 83 (col:Values must be of type <class 'str'>)
   C6: 18 (col:Values must be of type <class 'str'>)
   C7: '2; (Control=1, Dementia=2)' (col:Values must be of type <class 'str'>)
   C8: None (col:Values must be of type <class 'str'>)
   C10: 'X (24 subclass, 139 supertype)' (col:None)
   C11: 'X' (col:None)
Row 2:
   C2: 'X (10x v2 and v3)\nSample prep' (col:None)
   C3: 'Breast cancer [19]' (col:Values must be of type <class 'str'>)
   C4: '714 331' (col:Values must be of type <class 'str'>)
   C5: 126 (col:Values must be of type <class 'str'>)
   C6: 10 (col:Values must be of type <class 'str'>)
   C7: '2; (Control=1, Cancer=2)' (col:Values must be of type <class 'str'>)
   C8: None (col:Values must be of type <class 'str'>)
   C9: 'No cancer cells (because contralateral unaffected breast samples were taken)' (col:None)
   C10: 'X' (col:None)
   C11: 'X' (col:None)
Row 3:
   C2: 'X (collection site)' (col:None)
   C3: 'Covid-19 PBMC [5]' (col:Values must be of type <class 'str'>)
   C4: '993 171' (col:Values must be of type <class 'str'>)
   C5: 151 (col:Values must be of type <class 'str'>)
   C6: 10 (col:Values must be of type <class 'str'>)
   C7: '3; (Normal=1,Mild=2, Isch=3)' (col:Values must be of type <class 'str'>)
   C8: None (col:Values must be of type <class 'str'>)
   C10: 'X (64)' (col:None)
   C11: 'GEO' (col:None)
Row 4:
   C2: 'X (different studies)\n10x v2 and v3' (col:None)
   C3: 'Diabetes [28]' (col:FFFFFF00)
   C4: '264 235' (col:FFFFFF00)
   C5: 52 (col:FFFFFF00)
   C6: 13 (col:FFFFFF00)
   C7: '4; (Normal=1, EPD =2, Types I and II =3)' (col:FFFFFF00)
   C8: 'Possible but:\n- 9 dataset\n- Lots of conditions (Diabetes, chemical stress, strain, etc.)' (col:FFFFFF00)
   C10: 'X' (col:None)
   C11: 'X' (col:None)
Row 5:
   C1: 'X (subset kidney, scRNA-seq, by disease)' (col:None)
   C2: 'X (need to check, according to reviewer multiple sites)' (col:None)
   C3: 'Kidney [7]' (col:Values must be of type <class 'str'>)
   C4: '104 314' (col:Values must be of type <class 'str'>)
   C5: 45 (col:Values must be of type <class 'str'>)
   C6: 14 (col:Values must be of type <class 'str'>)
   C7: '3; (Normal=1,CKD=2, AKI=3)' (col:Values must be of type <class 'str'>)
   C8: None (col:Values must be of type <class 'str'>)
   C9: 'Subsetted to "Kidney" only' (col:None)
   C10: 'X' (col:None)
   C11: 'X' (col:None)
Row 6:
   C3: 'Kidney cancer [29]' (col:FFFFC000)
   C4: '50 236' (col:FFFFC000)
   C5: 17 (col:FFFFC000)
   C6: 6 (col:FFFFC000)
   C7: '2; (Normal=1, Cancer=2)' (col:FFFFC000)
   C8: 'Too few samples' (col:FFFFC000)
   C9: 'Separates by "tumor cells"\nAlso, very simple (probably all will score 1)' (col:FFFFC000)
   C10: None (col:Values must be of type <class 'str'>)
Row 7:
   C1: 'X' (col:None)
   C3: 'Lupus PBMC [30]' (col:Values must be of type <class 'str'>)
   C4: '1 263 676' (col:Values must be of type <class 'str'>)
   C5: 261 (col:Values must be of type <class 'str'>)
   C6: 11 (col:Values must be of type <class 'str'>)
   C7: '2; (Normal=1, SLE=2)' (col:Values must be of type <class 'str'>)
   C8: None (col:Values must be of type <class 'str'>)
   C9: 'Described cell type compositional AND gene expressional changes (IFN)' (col:None)
   C10: 'X' (col:None)
Row 8:
   C3: 'Lung [4]' (col:FFFFC000)
   C4: '941 504' (col:FFFFC000)
   C5: 165 (col:FFFFC000)
   C6: 12 (col:FFFFC000)
   C7: '5; (Normal=1, COPD=2, SCLC, LA, NSCLC=3)' (col:FFFFC000)
   C8: '"- Too many technical conditions (as described in paper)"' (col:FFFFC000)
   C10: None (col:Values must be of type <class 'str'>)
Row 9:
   C1: 'X' (col:None)
   C3: 'Myocardial infraction (2) [6]' (col:Values must be of type <class 'str'>)
   C4: '132 888' (col:Values must be of type <class 'str'>)
   C5: 23 (col:Values must be of type <class 'str'>)
   C6: 11 (col:Values must be of type <class 'str'>)
   C7: '3; (Normal=1, Isch=2, Fib=3)' (col:Values must be of type <class 'str'>)
   C8: 'Many bio conditions (location + disease)' (col:Values must be of type <class 'str'>)
   C9: None (col:Values must be of type <class 'str'>)
   C10: 'No (only clustering)' (col:None)
   C11: 'X' (col:None)
Row 10:
   C1: None (col:Values must be of type <class 'str'>)
   C2: None (col:Values must be of type <class 'str'>)
   C3: 'Myocardial infraction (1) [6]' (col:Values must be of type <class 'str'>)
   C4: '115 517' (col:Values must be of type <class 'str'>)
   C5: 20 (col:Values must be of type <class 'str'>)
   C6: 11 (col:Values must be of type <class 'str'>)
   C7: '2; (Normal=1, IZ=2)' (col:Values must be of type <class 'str'>)
   C8: None (col:Values must be of type <class 'str'>)
   C9: 'Both myocard are the same dataset' (col:Values must be of type <class 'str'>)
   C10: None (col:Values must be of type <class 'str'>)
   C11: 'X' (col:Values must be of type <class 'str'>)
Row 11:
   C3: 'Pancreas (PDAC) [31]' (col:FFFFC000)
   C4: '57 530' (col:FFFFC000)
   C5: 35 (col:FFFFC000)
   C6: 10 (col:FFFFC000)
   C7: '2; (Normal=1, Cancer=2)' (col:FFFFC000)
   C8: None (col:FFFFC000)
   C9: '- Separates by "Ductal 2" cells (which might be cancer cells)\n- Can be included but very simple (probably all will score 1)' (col:FFFFC000)
   C10: 'No' (col:None)
Row 12:
   C2: 'X' (col:None)
   C3: 'Parkinson [8]' (col:FFFFC000)
   C4: '2 096 155' (col:FFFFC000)
   C5: 97 (col:FFFFC000)
   C6: 11 (col:FFFFC000)
   C7: '2; (Normal=1, Diseased=2)' (col:FFFFC000)
   C8: '"- Only 12 cell types\n- 5 different brain regions"' (col:FFFFC000)
   C9: 'Cell type composition dominated by brain region\n- Could serve as negative control?' (col:None)
   C10: 'No' (col:None)
   C11: 'X' (col:None)


```

```python
import openpyxl

wb = openpyxl.load_workbook('ECODA_PAPER_DATASETS.xlsx', data_only=True)
ws = wb['Sheet1']

def get_color(cell):
    if cell.fill and cell.fill.fill_type:
        fg = cell.fill.fgColor
        if fg:
            if fg.type == 'rgb':
                return fg.rgb
            elif fg.type == 'theme':
                return f"theme:{fg.theme}"
            elif fg.type == 'indexed':
                return f"indexed:{fg.indexed}"
    return None

headers = [cell.value for cell in ws[1]]
print("Headers:", headers)

data = []
for row_idx in range(2, ws.max_row + 1):
    row = ws[row_idx]
    row_data = {}
    row_colors = {}
    for col_idx, cell in enumerate(row):
        h = headers[col_idx] if col_idx < len(headers) else f"Col_{col_idx+1}"
        row_data[h] = cell.value
        row_colors[h] = get_color(cell)
    data.append({'row': row_idx, 'values': row_data, 'colors': row_colors})

for item in data:
    print(f"\n--- Row {item['row']} ---")
    data_name = item['values'].get('Data')
    colors = set(c for c in item['colors'].values() if c is not None)
    print(f"Data Name: {data_name} | Colors in row: {colors}")
    for k, v in item['values'].items():
        if v is not None:
            c = item['colors'][k]
            print(f"  {k}: {repr(v)} (color: {c})")


```

```text
Headers: ['Benchmark', 'Batch effect', 'Data', 'No. cells', 'No. samples', 'No. cell_types', 'No. classes and Severity Scores', 'Possible exclusion criteria', 'Comment', 'Has high-granularity cell types', 'Available on cellxgene']

--- Row 2 ---
Data Name: Alzheimer [9] | Colors in row: {'theme:9'}
  Benchmark: 'X' (color: None)
  Data: 'Alzheimer [9]' (color: theme:9)
  No. cells: '1 395 601' (color: theme:9)
  No. samples: 83 (color: theme:9)
  No. cell_types: 18 (color: theme:9)
  No. classes and Severity Scores: '2; (Control=1, Dementia=2)' (color: theme:9)
  Has high-granularity cell types: 'X (24 subclass, 139 supertype)' (color: None)
  Available on cellxgene: 'X' (color: None)

--- Row 3 ---
Data Name: Breast cancer [19] | Colors in row: {'theme:9'}
  Batch effect: 'X (10x v2 and v3)\nSample prep' (color: None)
  Data: 'Breast cancer [19]' (color: theme:9)
  No. cells: '714 331' (color: theme:9)
  No. samples: 126 (color: theme:9)
  No. cell_types: 10 (color: theme:9)
  No. classes and Severity Scores: '2; (Control=1, Cancer=2)' (color: theme:9)
  Comment: 'No cancer cells (because contralateral unaffected breast samples were taken)' (color: None)
  Has high-granularity cell types: 'X' (color: None)
  Available on cellxgene: 'X' (color: None)

--- Row 4 ---
Data Name: Covid-19 PBMC [5] | Colors in row: {'theme:9'}
  Batch effect: 'X (collection site)' (color: None)
  Data: 'Covid-19 PBMC [5]' (color: theme:9)
  No. cells: '993 171' (color: theme:9)
  No. samples: 151 (color: theme:9)
  No. cell_types: 10 (color: theme:9)
  No. classes and Severity Scores: '3; (Normal=1,Mild=2, Isch=3)' (color: theme:9)
  Has high-granularity cell types: 'X (64)' (color: None)
  Available on cellxgene: 'GEO' (color: None)

--- Row 5 ---
Data Name: Diabetes [28] | Colors in row: {'FFFFFF00'}
  Batch effect: 'X (different studies)\n10x v2 and v3' (color: None)
  Data: 'Diabetes [28]' (color: FFFFFF00)
  No. cells: '264 235' (color: FFFFFF00)
  No. samples: 52 (color: FFFFFF00)
  No. cell_types: 13 (color: FFFFFF00)
  No. classes and Severity Scores: '4; (Normal=1, EPD =2, Types I and II =3)' (color: FFFFFF00)
  Possible exclusion criteria: 'Possible but:\n- 9 dataset\n- Lots of conditions (Diabetes, chemical stress, strain, etc.)' (color: FFFFFF00)
  Has high-granularity cell types: 'X' (color: None)
  Available on cellxgene: 'X' (color: None)

--- Row 6 ---
Data Name: Kidney [7] | Colors in row: {'theme:9'}
  Benchmark: 'X (subset kidney, scRNA-seq, by disease)' (color: None)
  Batch effect: 'X (need to check, according to reviewer multiple sites)' (color: None)
  Data: 'Kidney [7]' (color: theme:9)
  No. cells: '104 314' (color: theme:9)
  No. samples: 45 (color: theme:9)
  No. cell_types: 14 (color: theme:9)
  No. classes and Severity Scores: '3; (Normal=1,CKD=2, AKI=3)' (color: theme:9)
  Comment: 'Subsetted to "Kidney" only' (color: None)
  Has high-granularity cell types: 'X' (color: None)
  Available on cellxgene: 'X' (color: None)

--- Row 7 ---
Data Name: Kidney cancer [29] | Colors in row: {'theme:1', 'FFFFC000'}
  Data: 'Kidney cancer [29]' (color: FFFFC000)
  No. cells: '50 236' (color: FFFFC000)
  No. samples: 17 (color: FFFFC000)
  No. cell_types: 6 (color: FFFFC000)
  No. classes and Severity Scores: '2; (Normal=1, Cancer=2)' (color: FFFFC000)
  Possible exclusion criteria: 'Too few samples' (color: FFFFC000)
  Comment: 'Separates by "tumor cells"\nAlso, very simple (probably all will score 1)' (color: FFFFC000)

--- Row 8 ---
Data Name: Lupus PBMC [30] | Colors in row: {'theme:9'}
  Benchmark: 'X' (color: None)
  Data: 'Lupus PBMC [30]' (color: theme:9)
  No. cells: '1 263 676' (color: theme:9)
  No. samples: 261 (color: theme:9)
  No. cell_types: 11 (color: theme:9)
  No. classes and Severity Scores: '2; (Normal=1, SLE=2)' (color: theme:9)
  Comment: 'Described cell type compositional AND gene expressional changes (IFN)' (color: None)
  Has high-granularity cell types: 'X' (color: None)

--- Row 9 ---
Data Name: Lung [4] | Colors in row: {'theme:1', 'FFFFC000'}
  Data: 'Lung [4]' (color: FFFFC000)
  No. cells: '941 504' (color: FFFFC000)
  No. samples: 165 (color: FFFFC000)
  No. cell_types: 12 (color: FFFFC000)
  No. classes and Severity Scores: '5; (Normal=1, COPD=2, SCLC, LA, NSCLC=3)' (color: FFFFC000)
  Possible exclusion criteria: '"- Too many technical conditions (as described in paper)"' (color: FFFFC000)

--- Row 10 ---
Data Name: Myocardial infraction (2) [6] | Colors in row: {'theme:9'}
  Benchmark: 'X' (color: None)
  Data: 'Myocardial infraction (2) [6]' (color: theme:9)
  No. cells: '132 888' (color: theme:9)
  No. samples: 23 (color: theme:9)
  No. cell_types: 11 (color: theme:9)
  No. classes and Severity Scores: '3; (Normal=1, Isch=2, Fib=3)' (color: theme:9)
  Possible exclusion criteria: 'Many bio conditions (location + disease)' (color: theme:9)
  Has high-granularity cell types: 'No (only clustering)' (color: None)
  Available on cellxgene: 'X' (color: None)

--- Row 11 ---
Data Name: Myocardial infraction (1) [6] | Colors in row: {'theme:0'}
  Data: 'Myocardial infraction (1) [6]' (color: theme:0)
  No. cells: '115 517' (color: theme:0)
  No. samples: 20 (color: theme:0)
  No. cell_types: 11 (color: theme:0)
  No. classes and Severity Scores: '2; (Normal=1, IZ=2)' (color: theme:0)
  Comment: 'Both myocard are the same dataset' (color: theme:0)
  Available on cellxgene: 'X' (color: theme:0)

--- Row 12 ---
Data Name: Pancreas (PDAC) [31] | Colors in row: {'FFFFC000'}
  Data: 'Pancreas (PDAC) [31]' (color: FFFFC000)
  No. cells: '57 530' (color: FFFFC000)
  No. samples: 35 (color: FFFFC000)
  No. cell_types: 10 (color: FFFFC000)
  No. classes and Severity Scores: '2; (Normal=1, Cancer=2)' (color: FFFFC000)
  Comment: '- Separates by "Ductal 2" cells (which might be cancer cells)\n- Can be included but very simple (probably all will score 1)' (color: FFFFC000)
  Has high-granularity cell types: 'No' (color: None)

--- Row 13 ---
Data Name: Parkinson [8] | Colors in row: {'FFFFC000'}
  Batch effect: 'X' (color: None)
  Data: 'Parkinson [8]' (color: FFFFC000)
  No. cells: '2 096 155' (color: FFFFC000)
  No. samples: 97 (color: FFFFC000)
  No. cell_types: 11 (color: FFFFC000)
  No. classes and Severity Scores: '2; (Normal=1, Diseased=2)' (color: FFFFC000)
  Possible exclusion criteria: '"- Only 12 cell types\n- 5 different brain regions"' (color: FFFFC000)
  Comment: 'Cell type composition dominated by brain region\n- Could serve as negative control?' (color: None)
  Has high-granularity cell types: 'No' (color: None)
  Available on cellxgene: 'X' (color: None)


```

```python
# Let's inspect the exact color hex or theme color mapping for each row

import openpyxl

wb = openpyxl.load_workbook('ECODA_PAPER_DATASETS.xlsx', data_only=True)
ws = wb['Sheet1']

for r in range(2, ws.max_row + 1):
    c_data = ws.cell(row=r, column=3)
    fill = c_data.fill
    fg = fill.fgColor
    print(f"Row {r:2d} | Dataset: {c_data.value:<30} | fill_type: {fill.fill_type} | fgColor type: {fg.type} | rgb: {fg.rgb} | theme: {fg.theme} | tint: {fg.tint}")


```

```text
Row  2 | Dataset: Alzheimer [9]                  | fill_type: solid | fgColor type: theme | rgb: Values must be of type <class 'str'> | theme: 9 | tint: 0.3999755851924192
Row  3 | Dataset: Breast cancer [19]             | fill_type: solid | fgColor type: theme | rgb: Values must be of type <class 'str'> | theme: 9 | tint: 0.7999816888943144
Row  4 | Dataset: Covid-19 PBMC [5]              | fill_type: solid | fgColor type: theme | rgb: Values must be of type <class 'str'> | theme: 9 | tint: 0.7999816888943144
Row  5 | Dataset: Diabetes [28]                  | fill_type: solid | fgColor type: rgb | rgb: FFFFFF00 | theme: Values must be of type <class 'int'> | tint: 0.0
Row  6 | Dataset: Kidney [7]                     | fill_type: solid | fgColor type: theme | rgb: Values must be of type <class 'str'> | theme: 9 | tint: 0.3999755851924192
Row  7 | Dataset: Kidney cancer [29]             | fill_type: solid | fgColor type: rgb | rgb: FFFFC000 | theme: Values must be of type <class 'int'> | tint: 0.0
Row  8 | Dataset: Lupus PBMC [30]                | fill_type: solid | fgColor type: theme | rgb: Values must be of type <class 'str'> | theme: 9 | tint: 0.3999755851924192
Row  9 | Dataset: Lung [4]                       | fill_type: solid | fgColor type: rgb | rgb: FFFFC000 | theme: Values must be of type <class 'int'> | tint: 0.0
Row 10 | Dataset: Myocardial infraction (2) [6]  | fill_type: solid | fgColor type: theme | rgb: Values must be of type <class 'str'> | theme: 9 | tint: 0.7999816888943144
Row 11 | Dataset: Myocardial infraction (1) [6]  | fill_type: solid | fgColor type: theme | rgb: Values must be of type <class 'str'> | theme: 0 | tint: -0.249977111117893
Row 12 | Dataset: Pancreas (PDAC) [31]           | fill_type: solid | fgColor type: rgb | rgb: FFFFC000 | theme: Values must be of type <class 'int'> | tint: 0.0
Row 13 | Dataset: Parkinson [8]                  | fill_type: solid | fgColor type: rgb | rgb: FFFFC000 | theme: Values must be of type <class 'int'> | tint: 0.0


```

```python
import json
import re

# Map theme/rgb to feasibility status
# theme 9 = Green
# FFFFFF00 = Yellow
# FFFFC000 = Orange
# theme 0 = Duplicate / Uncolored (Myocardial infraction (1) [6])

datasets = []

def parse_num_cells(val):
    if val is None:
        return None
    val_str = str(val).replace(' ', '').strip()
    try:
        return int(val_str)
    except:
        return val_str

def parse_ref_and_name(val):
    if not val:
        return None, None
    m = re.search(r'^(.*?)\s*\[(\d+)\]$', str(val).strip())
    if m:
        return m.group(1).strip(), int(m.group(2))
    return str(val).strip(), None

for r in range(2, ws.max_row + 1):
    c_data = ws.cell(row=r, column=3)
    data_raw = c_data.value
    name, ref_id = parse_ref_and_name(data_raw)
    
    fill = c_data.fill
    fg = fill.fgColor
    
    if fg.type == 'rgb':
        color_code = fg.rgb
    elif fg.type == 'theme':
        color_code = f"theme_{fg.theme}"
    else:
        color_code = "none"
        
    if color_code == 'theme_9':
        status = "Green (Good, feasible)"
        feasibility_code = "green"
    elif color_code == 'FFFFFF00':
        status = "Yellow (Might need to be checked)"
        feasibility_code = "yellow"
    elif color_code == 'FFFFC000':
        status = "Orange (Probably not, needs to be checked)"
        feasibility_code = "orange"
    else:
        status = "Neutral / Duplicate"
        feasibility_code = "neutral"

    row_dict = {
        "dataset_name": name,
        "paper_reference_id": ref_id,
        "feasibility_status": status,
        "feasibility_code": feasibility_code,
        "benchmark": ws.cell(row=r, column=1).value,
        "batch_effect": ws.cell(row=r, column=2).value,
        "no_cells": parse_num_cells(ws.cell(row=r, column=4).value),
        "no_samples": ws.cell(row=r, column=5).value,
        "no_cell_types": ws.cell(row=r, column=6).value,
        "classes_and_severity_scores": ws.cell(row=r, column=7).value,
        "exclusion_criteria": ws.cell(row=r, column=8).value,
        "comment": ws.cell(row=r, column=9).value,
        "has_high_granularity_cell_types": ws.cell(row=r, column=10).value,
        "available_on_cellxgene": ws.cell(row=r, column=11).value
    }
    datasets.append(row_dict)

print(json.dumps(datasets, indent=2))


```

```text
[
  {
    "dataset_name": "Alzheimer",
    "paper_reference_id": 9,
    "feasibility_status": "Green (Good, feasible)",
    "feasibility_code": "green",
    "benchmark": "X",
    "batch_effect": null,
    "no_cells": 1395601,
    "no_samples": 83,
    "no_cell_types": 18,
    "classes_and_severity_scores": "2; (Control=1, Dementia=2)",
    "exclusion_criteria": null,
    "comment": null,
    "has_high_granularity_cell_types": "X (24 subclass, 139 supertype)",
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Breast cancer",
    "paper_reference_id": 19,
    "feasibility_status": "Green (Good, feasible)",
    "feasibility_code": "green",
    "benchmark": null,
    "batch_effect": "X (10x v2 and v3)\nSample prep",
    "no_cells": 714331,
    "no_samples": 126,
    "no_cell_types": 10,
    "classes_and_severity_scores": "2; (Control=1, Cancer=2)",
    "exclusion_criteria": null,
    "comment": "No cancer cells (because contralateral unaffected breast samples were taken)",
    "has_high_granularity_cell_types": "X",
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Covid-19 PBMC",
    "paper_reference_id": 5,
    "feasibility_status": "Green (Good, feasible)",
    "feasibility_code": "green",
    "benchmark": null,
    "batch_effect": "X (collection site)",
    "no_cells": 993171,
    "no_samples": 151,
    "no_cell_types": 10,
    "classes_and_severity_scores": "3; (Normal=1,Mild=2, Isch=3)",
    "exclusion_criteria": null,
    "comment": null,
    "has_high_granularity_cell_types": "X (64)",
    "available_on_cellxgene": "GEO"
  },
  {
    "dataset_name": "Diabetes",
    "paper_reference_id": 28,
    "feasibility_status": "Yellow (Might need to be checked)",
    "feasibility_code": "yellow",
    "benchmark": null,
    "batch_effect": "X (different studies)\n10x v2 and v3",
    "no_cells": 264235,
    "no_samples": 52,
    "no_cell_types": 13,
    "classes_and_severity_scores": "4; (Normal=1, EPD =2, Types I and II =3)",
    "exclusion_criteria": "Possible but:\n- 9 dataset\n- Lots of conditions (Diabetes, chemical stress, strain, etc.)",
    "comment": null,
    "has_high_granularity_cell_types": "X",
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Kidney",
    "paper_reference_id": 7,
    "feasibility_status": "Green (Good, feasible)",
    "feasibility_code": "green",
    "benchmark": "X (subset kidney, scRNA-seq, by disease)",
    "batch_effect": "X (need to check, according to reviewer multiple sites)",
    "no_cells": 104314,
    "no_samples": 45,
    "no_cell_types": 14,
    "classes_and_severity_scores": "3; (Normal=1,CKD=2, AKI=3)",
    "exclusion_criteria": null,
    "comment": "Subsetted to \"Kidney\" only",
    "has_high_granularity_cell_types": "X",
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Kidney cancer",
    "paper_reference_id": 29,
    "feasibility_status": "Orange (Probably not, needs to be checked)",
    "feasibility_code": "orange",
    "benchmark": null,
    "batch_effect": null,
    "no_cells": 50236,
    "no_samples": 17,
    "no_cell_types": 6,
    "classes_and_severity_scores": "2; (Normal=1, Cancer=2)",
    "exclusion_criteria": "Too few samples",
    "comment": "Separates by \"tumor cells\"\nAlso, very simple (probably all will score 1)",
    "has_high_granularity_cell_types": null,
    "available_on_cellxgene": null
  },
  {
    "dataset_name": "Lupus PBMC",
    "paper_reference_id": 30,
    "feasibility_status": "Green (Good, feasible)",
    "feasibility_code": "green",
    "benchmark": "X",
    "batch_effect": null,
    "no_cells": 1263676,
    "no_samples": 261,
    "no_cell_types": 11,
    "classes_and_severity_scores": "2; (Normal=1, SLE=2)",
    "exclusion_criteria": null,
    "comment": "Described cell type compositional AND gene expressional changes (IFN)",
    "has_high_granularity_cell_types": "X",
    "available_on_cellxgene": null
  },
  {
    "dataset_name": "Lung",
    "paper_reference_id": 4,
    "feasibility_status": "Orange (Probably not, needs to be checked)",
    "feasibility_code": "orange",
    "benchmark": null,
    "batch_effect": null,
    "no_cells": 941504,
    "no_samples": 165,
    "no_cell_types": 12,
    "classes_and_severity_scores": "5; (Normal=1, COPD=2, SCLC, LA, NSCLC=3)",
    "exclusion_criteria": "\"- Too many technical conditions (as described in paper)\"",
    "comment": null,
    "has_high_granularity_cell_types": null,
    "available_on_cellxgene": null
  },
  {
    "dataset_name": "Myocardial infraction (2)",
    "paper_reference_id": 6,
    "feasibility_status": "Green (Good, feasible)",
    "feasibility_code": "green",
    "benchmark": "X",
    "batch_effect": null,
    "no_cells": 132888,
    "no_samples": 23,
    "no_cell_types": 11,
    "classes_and_severity_scores": "3; (Normal=1, Isch=2, Fib=3)",
    "exclusion_criteria": "Many bio conditions (location + disease)",
    "comment": null,
    "has_high_granularity_cell_types": "No (only clustering)",
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Myocardial infraction (1)",
    "paper_reference_id": 6,
    "feasibility_status": "Neutral / Duplicate",
    "feasibility_code": "neutral",
    "benchmark": null,
    "batch_effect": null,
    "no_cells": 115517,
    "no_samples": 20,
    "no_cell_types": 11,
    "classes_and_severity_scores": "2; (Normal=1, IZ=2)",
    "exclusion_criteria": null,
    "comment": "Both myocard are the same dataset",
    "has_high_granularity_cell_types": null,
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Pancreas (PDAC)",
    "paper_reference_id": 31,
    "feasibility_status": "Orange (Probably not, needs to be checked)",
    "feasibility_code": "orange",
    "benchmark": null,
    "batch_effect": null,
    "no_cells": 57530,
    "no_samples": 35,
    "no_cell_types": 10,
    "classes_and_severity_scores": "2; (Normal=1, Cancer=2)",
    "exclusion_criteria": null,
    "comment": "- Separates by \"Ductal 2\" cells (which might be cancer cells)\n- Can be included but very simple (probably all will score 1)",
    "has_high_granularity_cell_types": "No",
    "available_on_cellxgene": null
  },
  {
    "dataset_name": "Parkinson",
    "paper_reference_id": 8,
    "feasibility_status": "Orange (Probably not, needs to be checked)",
    "feasibility_code": "orange",
    "benchmark": null,
    "batch_effect": "X",
    "no_cells": 2096155,
    "no_samples": 97,
    "no_cell_types": 11,
    "classes_and_severity_scores": "2; (Normal=1, Diseased=2)",
    "exclusion_criteria": "\"- Only 12 cell types\n- 5 different brain regions\"",
    "comment": "Cell type composition dominated by brain region\n- Could serve as negative control?",
    "has_high_granularity_cell_types": "No",
    "available_on_cellxgene": "X"
  }
]


```

Here is the converted dataset in **machine/agent-readable JSON format**, along with an **explicit Markdown summary table**.

The color coding from the Excel spreadsheet has been extracted into explicit status fields (`feasibility_status` and `feasibility_code`):

* **`green`**: Good, should be feasible
* **`yellow`**: Might need to be checked
* **`orange`**: Probably not, but needs to be checked
* **`neutral`**: Duplicate/related note entry

---

### 1. Agent-Readable JSON Format

```json
[
  {
    "dataset_name": "Alzheimer",
    "paper_reference_id": 9,
    "feasibility_code": "green",
    "feasibility_status": "Good, should be feasible",
    "benchmark": "X",
    "batch_effect": null,
    "no_cells": 1395601,
    "no_samples": 83,
    "no_cell_types": 18,
    "classes_and_severity_scores": "2; (Control=1, Dementia=2)",
    "exclusion_criteria": null,
    "comment": null,
    "has_high_granularity_cell_types": "X (24 subclass, 139 supertype)",
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Breast cancer",
    "paper_reference_id": 19,
    "feasibility_code": "green",
    "feasibility_status": "Good, should be feasible",
    "benchmark": null,
    "batch_effect": "X (10x v2 and v3)\nSample prep",
    "no_cells": 714331,
    "no_samples": 126,
    "no_cell_types": 10,
    "classes_and_severity_scores": "2; (Control=1, Cancer=2)",
    "exclusion_criteria": null,
    "comment": "No cancer cells (because contralateral unaffected breast samples were taken)",
    "has_high_granularity_cell_types": "X",
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Covid-19 PBMC",
    "paper_reference_id": 5,
    "feasibility_code": "green",
    "feasibility_status": "Good, should be feasible",
    "benchmark": null,
    "batch_effect": "X (collection site)",
    "no_cells": 993171,
    "no_samples": 151,
    "no_cell_types": 10,
    "classes_and_severity_scores": "3; (Normal=1,Mild=2, Isch=3)",
    "exclusion_criteria": null,
    "comment": null,
    "has_high_granularity_cell_types": "X (64)",
    "available_on_cellxgene": "GEO"
  },
  {
    "dataset_name": "Diabetes",
    "paper_reference_id": 28,
    "feasibility_code": "yellow",
    "feasibility_status": "Might need to be checked",
    "benchmark": null,
    "batch_effect": "X (different studies)\n10x v2 and v3",
    "no_cells": 264235,
    "no_samples": 52,
    "no_cell_types": 13,
    "classes_and_severity_scores": "4; (Normal=1, EPD =2, Types I and II =3)",
    "exclusion_criteria": "Possible but:\n- 9 dataset\n- Lots of conditions (Diabetes, chemical stress, strain, etc.)",
    "comment": null,
    "has_high_granularity_cell_types": "X",
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Kidney",
    "paper_reference_id": 7,
    "feasibility_code": "green",
    "feasibility_status": "Good, should be feasible",
    "benchmark": "X (subset kidney, scRNA-seq, by disease)",
    "batch_effect": "X (need to check, according to reviewer multiple sites)",
    "no_cells": 104314,
    "no_samples": 45,
    "no_cell_types": 14,
    "classes_and_severity_scores": "3; (Normal=1,CKD=2, AKI=3)",
    "exclusion_criteria": null,
    "comment": "Subsetted to \"Kidney\" only",
    "has_high_granularity_cell_types": "X",
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Kidney cancer",
    "paper_reference_id": 29,
    "feasibility_code": "orange",
    "feasibility_status": "Probably not, needs to be checked",
    "benchmark": null,
    "batch_effect": null,
    "no_cells": 50236,
    "no_samples": 17,
    "no_cell_types": 6,
    "classes_and_severity_scores": "2; (Normal=1, Cancer=2)",
    "exclusion_criteria": "Too few samples",
    "comment": "Separates by \"tumor cells\"\nAlso, very simple (probably all will score 1)",
    "has_high_granularity_cell_types": null,
    "available_on_cellxgene": null
  },
  {
    "dataset_name": "Lupus PBMC",
    "paper_reference_id": 30,
    "feasibility_code": "green",
    "feasibility_status": "Good, should be feasible",
    "benchmark": "X",
    "batch_effect": null,
    "no_cells": 1263676,
    "no_samples": 261,
    "no_cell_types": 11,
    "classes_and_severity_scores": "2; (Normal=1, SLE=2)",
    "exclusion_criteria": null,
    "comment": "Described cell type compositional AND gene expressional changes (IFN)",
    "has_high_granularity_cell_types": "X",
    "available_on_cellxgene": null
  },
  {
    "dataset_name": "Lung",
    "paper_reference_id": 4,
    "feasibility_code": "orange",
    "feasibility_status": "Probably not, needs to be checked",
    "benchmark": null,
    "batch_effect": null,
    "no_cells": 941504,
    "no_samples": 165,
    "no_cell_types": 12,
    "classes_and_severity_scores": "5; (Normal=1, COPD=2, SCLC, LA, NSCLC=3)",
    "exclusion_criteria": "- Too many technical conditions (as described in paper)",
    "comment": null,
    "has_high_granularity_cell_types": null,
    "available_on_cellxgene": null
  },
  {
    "dataset_name": "Myocardial infraction (2)",
    "paper_reference_id": 6,
    "feasibility_code": "green",
    "feasibility_status": "Good, should be feasible",
    "benchmark": "X",
    "batch_effect": null,
    "no_cells": 132888,
    "no_samples": 23,
    "no_cell_types": 11,
    "classes_and_severity_scores": "3; (Normal=1, Isch=2, Fib=3)",
    "exclusion_criteria": "Many bio conditions (location + disease)",
    "comment": null,
    "has_high_granularity_cell_types": "No (only clustering)",
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Myocardial infraction (1)",
    "paper_reference_id": 6,
    "feasibility_code": "neutral",
    "feasibility_status": "Neutral / Duplicate entry note",
    "benchmark": null,
    "batch_effect": null,
    "no_cells": 115517,
    "no_samples": 20,
    "no_cell_types": 11,
    "classes_and_severity_scores": "2; (Normal=1, IZ=2)",
    "exclusion_criteria": null,
    "comment": "Both myocard are the same dataset",
    "has_high_granularity_cell_types": null,
    "available_on_cellxgene": "X"
  },
  {
    "dataset_name": "Pancreas (PDAC)",
    "paper_reference_id": 31,
    "feasibility_code": "orange",
    "feasibility_status": "Probably not, needs to be checked",
    "benchmark": null,
    "batch_effect": null,
    "no_cells": 57530,
    "no_samples": 35,
    "no_cell_types": 10,
    "classes_and_severity_scores": "2; (Normal=1, Cancer=2)",
    "exclusion_criteria": null,
    "comment": "- Separates by \"Ductal 2\" cells (which might be cancer cells)\n- Can be included but very simple (probably all will score 1)",
    "has_high_granularity_cell_types": "No",
    "available_on_cellxgene": null
  },
  {
    "dataset_name": "Parkinson",
    "paper_reference_id": 8,
    "feasibility_code": "orange",
    "feasibility_status": "Probably not, needs to be checked",
    "benchmark": null,
    "batch_effect": "X",
    "no_cells": 2096155,
    "no_samples": 97,
    "no_cell_types": 11,
    "classes_and_severity_scores": "2; (Normal=1, Diseased=2)",
    "exclusion_criteria": "- Only 12 cell types\n- 5 different brain regions",
    "comment": "Cell type composition dominated by brain region\n- Could serve as negative control?",
    "has_high_granularity_cell_types": "No",
    "available_on_cellxgene": "X"
  }
]

```

---

### 2. Summary Table by Feasibility Category

| Feasibility | Dataset Name | Ref ID | No. Cells | No. Samples | Classes & Severity | Main Issues / Exclusion Criteria |
| --- | --- | --- | --- | --- | --- | --- |
| **🟢 Green** (Good) | **Alzheimer** | [9] | 1,395,601 | 83 | 2 (Control=1, Dementia=2) | High-granularity (24 subclass, 139 supertype) |
| **🟢 Green** (Good) | **Breast cancer** | [19] | 714,331 | 126 | 2 (Control=1, Cancer=2) | No cancer cells (unaffected breast samples) |
| **🟢 Green** (Good) | **Covid-19 PBMC** | [5] | 993,171 | 151 | 3 (Normal=1, Mild=2, Isch=3) | Available on GEO |
| **🟢 Green** (Good) | **Kidney** | [7] | 104,314 | 45 | 3 (Normal=1, CKD=2, AKI=3) | Subsetted to Kidney only; multi-site batch effect |
| **🟢 Green** (Good) | **Lupus PBMC** | [30] | 1,263,676 | 261 | 2 (Normal=1, SLE=2) | Cell type composition & IFN gene expression |
| **🟢 Green** (Good) | **Myocardial infraction (2)** | [6] | 132,888 | 23 | 3 (Normal=1, Isch=2, Fib=3) | Many bio conditions (location + disease) |
| **🟡 Yellow** (Check) | **Diabetes** | [28] | 264,235 | 52 | 4 (Normal=1, EPD=2, Types I/II=3) | 9 sub-datasets; many technical/bio conditions |
| **🟠 Orange** (Unlikely) | **Kidney cancer** | [29] | 50,236 | 17 | 2 (Normal=1, Cancer=2) | Too few samples; very simple separation |
| **🟠 Orange** (Unlikely) | **Lung** | [4] | 941,504 | 165 | 5 (Normal=1, COPD=2, SCLC, LA, NSCLC=3) | Too many technical conditions |
| **🟠 Orange** (Unlikely) | **Pancreas (PDAC)** | [31] | 57,530 | 35 | 2 (Normal=1, Cancer=2) | Separates by Ductal 2 cells; very simple |
| **🟠 Orange** (Unlikely) | **Parkinson** | [8] | 2,096,155 | 97 | 2 (Normal=1, Diseased=2) | 5 brain regions dominate composition |
| **⚪ Neutral** (Dup) | **Myocardial infraction (1)** | [6] | 115,517 | 20 | 2 (Normal=1, IZ=2) | Duplicate dataset of Myocardial Infarction (2) |