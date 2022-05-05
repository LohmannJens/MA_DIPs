'''
Loads the data for the deletion sides (start/end point) from Alnaji 2019.
Takes a look at the nucleotide distribution around these points.
Goal is to see if any of the four nucleotides occure more/less often.
'''
import os



from utils import get_seq_len, load_excel

file_path = os.path.join("data", "DI_Influenza_FA_JVI.xlsx")
cleaned_data_dict = load_excel(file_path)


print(cleaned_data_dict)
