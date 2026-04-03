import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from itertools import chain

tree = ET.parse('data/icd10cm-tabular-April-2024.xml')
root = tree.getroot()
diags = root.findall(".//diag")
table_icd10 = np.array([[diag.findall("name")[0].text, diag.findall("desc")[0].text] for diag in diags])
table_icd10 = pd.DataFrame(table_icd10, columns=["code", "description"]).set_index("code")

# icd10.table_icd10