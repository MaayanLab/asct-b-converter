#!/usr/bin/env python
# coding: utf-8

# ## Harmonize ETL: HuBMAP
# 
# Data Source: https://hubmapconsortium.github.io/ccf-asct-reporter/vis?selectedOrgans=All_Organs-v1.1,Blood-v1.1,Blood_Vasculature-v1.1,bone_marrow-v1.1,brain-v1.1,Eye-v1.0,Fallopian_Tube-v1.1,heart-v1.1,kidney-v1.1,Knee-v1.0,large_intestine-v1.1,Liver-v1.0,lung-v1.1,lymph_nodes-v1.1,Lymph_Vasculature-v1.0,Ovary-v1.0,Pancreas-v1.0,Peripheral_Nervous_System-v1.0,Prostate-v1.0,skin-v1.1,Small_Intestine-v1.0,spleen-v1.1,thymus-v1.1,Ureter-v1.0,Urinary_Bladder-v1.0,Uterus-v1.0&playground=false

# In[1]:


import re
import sys
import pprint 
from datetime import date


# ### Notebook Information

# In[2]:


print('This notebook was run on:', date.today(), '\n Python version:', sys.version)


# # Initialization

# In[3]:


# Enter input file name
file_name = "Uterus.tsv"


# In[4]:


read_file = open(file_name, "r")


# In[5]:


def getList(dict):
    return dict.keys()


# ### Create a dictionary where the keys are columns (CT/BGene) and values are rows

# In[6]:


# Create dictionary where keys are columns and values are rows

def tsv_to_columns_row_dictionary(read_file):
    excel_dict = dict()
    key_list = list()
    for count, row in enumerate(read_file):
        temp = re.split(r'\t', row.strip("\n"))
        if count == 8:
            version = temp[1]
        if count == 10:
            for value in temp:
                excel_dict[value] = list()
                key_list.append(value)       
        if count > 10:
            for count, value in enumerate(temp):
                excel_dict[key_list[count]].append(value)
    return(excel_dict, key_list, version)


# ## Output of column/row dictionary

# In[7]:


excel_dict, key_list, version = tsv_to_columns_row_dictionary(read_file)


# ### Create a list of cell types. Cell types are listed at the highest resolution base on ASCT+B Table

# In[8]:


# Key_list is a list of column name. Ex. CT/2... 
# Add cell type to dictionary as key
def add_cell_type_to_list(excel_dict, key_list):
    cell_type_list = list()
    for index in range(0, len(excel_dict["AS/1"])):
        cell_type_list.append(list())
    
    for index in range(0, len(excel_dict["AS/1"])):
        for label in key_list:
            if "CT" in label and "LABEL" not in label and "ID" not in label and "NOTES" not in label:
                if excel_dict[label][index]:
                    cell_type_list[index] = (excel_dict[label][index])
    
    # Add unique identifier
    final_cell_type_list = list()
    for count, cell_type in enumerate(cell_type_list):
        final_cell_type = "{} ({})".format(cell_type, count)
        final_cell_type_list.append(final_cell_type)
    
    
    return(final_cell_type_list)


# ### Cell Types List

# In[9]:


cell_type_list = add_cell_type_to_list(excel_dict, key_list)
#pprint.pprint(cell_type_list)


# ### Create a dictionary where the cell type is the key and the value is a gene marker

# In[10]:


def add_gene_marker(cell_type_list, excel_dict):
    cell_type_gene_marker_dict = dict()
    
    for count, cell_type in enumerate(cell_type_list):
        cell_type_gene_marker_dict[cell_type] = excel_dict["All Gene Biomarkers"][count]
    
    return(cell_type_gene_marker_dict)


def add_gene_marker_no_all(cell_type_list, excel_dict):
    gene_list = list()
    for index in range(0, len(excel_dict["AS/1"])):
        gene_list.append(list())
    
    for index in range(0, len(excel_dict["AS/1"])):
        for label in key_list:
            if "BGene" in label and "LABEL" not in label and "ID" not in label:
                if excel_dict[label][index]:
                    if gene_list[index]:                       
                        gene_list[index].append(excel_dict[label][index])
                    else:
                        gene_list[index] = list()
                        gene_list[index].append(excel_dict[label][index])
                    
    cell_type_gene_marker_dict = dict()
    
    for count, cell_type in enumerate(cell_type_list):
        cell_type_gene_marker_dict[cell_type] = gene_list[count]
    
    return(cell_type_gene_marker_dict)


# ## Output of cell type/gene marker dictionary

# In[11]:


try:
    final_dict = add_gene_marker(cell_type_list, excel_dict)
except:
    final_dict = add_gene_marker_no_all(cell_type_list, excel_dict)
pprint.pprint(final_dict)


# In[12]:


output_filename = "{}_GMT_{}.tsv".format(file_name[:-4], version)
fh = open(output_filename , "w")
for key, value in final_dict.items():
    if "[" in key and not value:
        continue
    elif value and type(value) == list:
        fh.write("{}\t{}\n".format(key.strip('\t'), ', '.join(value).strip('\t')))
    elif value:
        fh.write("{}\t{}\n".format(key, value))
    else:
        fh.write("{}\n".format(key))


# ### Generate url to download ASCT+B Table

# In[13]:


sheetID_dict = {"Blood": "1ZYcSWnFHmzR9XKy_002f_oA4PfzokiW4IxkaZZOusvg",
                "Blood_Vascular": "1IlELzPwpWoHUcDAmNBWofXfislAaF_oR8yVpwy-zl18",
                "Bone_Marrow": "1tnqtCAWSA6atiUBUOOjAHdOrjDw_fsIoCd5RkAmw310",
                "Brain": "1TiwW1NZJ5kdCzJ4zwCpY3Gzv3WE5WUoBDWIAkU5gXd0",
                "Breast": "1Ac7C4dX7eYSMyR75AA2uVY9ZgNGOZZgbqgR8wmp-wdk",
                "Eye": "1u7IbxnPABRpYL5rFxOba8cmlvG1yGp-dwD3TV3V26K4",
                "Fallopian_Tube": "16tAvAmjwKwbq5SDz7UZ-T1N_KUHRGqPDbMqffFuInMI",
                "Heart": "1UhEZpDxQLCJLLx0gnWYDMQP8M-dwjZo_vIyPfjBCcVM",
                "Kidney": "1PgjYp4MEWANfbxGIxFsJ9vkfEU90MP-v3p5oVlH8U-E",
                "Knee": "1QidDho8DxBYjsxaqApiIZA__Z7aWnB61KvC422g2kx8",
                "Large_Intestine": "1vU6mQmnzAAxctbNYPoFxJ8NvbUql8pbipsGdt7YCOQQ",
                "Liver": "1tPDKw_znxqWhZYPTeVN4AN2_F4-JecsdeUgp2lj4P8g",
                "Lung": "1tK916JyG5ZSXW_cXfsyZnzXfjyoN-8B2GXLbYD6_vF0",
                "Lymph_Node": "1aK9gJ2_kMb2B8zrQgScDgxpEWAcCs7kl-gnQGwV3LHM",
                "Lymph_Vasculature": "1SILRNUI71BEVWl1fpsi_32DSuSA-bAPgXv5pTfKnrOE",
                "Ovary": "1FE2XufrruExUWqcai3XRFqtMjeEdzoLKJ-YNa-nRZ1M",
                "Pancreas": "1CIWqIygz2OzxMECIvhudFN14Kt7-JFUBLpzn5uuH5Xs",
                "Peripheral_Nervous_System": "1KifiEDn3PpJ8pjz9_ka4TWkT085wLIzIQP5NKSvb2Ac",
                "Placenta": "1TqatRIsZZ5QwvWdz6H4Un-sukbzSd21_x41Gqnn5UEY",
                "Prostate": "1_O5yXOesG93dobMHRSIvVAt9xj7mDnEAYdRJcHYJ84U",
                "Skin": "1Pmi3g26vhbg9HU6GDpIvxKbIP985JM-5eytOHxJUdZs",
                "Small_Intestine": "1Xlds8FzZ8ecmy3cxYJt1ijQC9FifamZRZ5KzH4Yt-MQ",
                "Spleen": "1HL7aHx5A2KOa1KsJ0PIagqxdshVavFIEJZP6_YDtUww",
                "Spinal_Cord": "1tK916JyG5ZSXW_cXfsyZnzXfjyoN-8B2GXLbYD6_vF0",
                "Thymus": "1nSiz2yFDMJSqIXbnAP_EXIQZfN6ZflOs-WBdZ6LVhUY",
                "Ureter": "1tK916JyG5ZSXW_cXfsyZnzXfjyoN-8B2GXLbYD6_vF0",
                "Urinary_Bladder": "1ohOG5jMf9d9eqjbVK6_u3CvgfG3wcLfs_pxB2838wOo",
                "Uterus": "1yEcbJMrUIzJY-4JNtF1Y_eUpAQsgKF6DX2-5Z3UXBeE"
               }


# In[14]:


gidID_dict = {"Blood": "360436225",
              "Blood_Vascular": "997949803",
              "Bone_Marrow": "771476671",
              "Brain": "2056967441",
              "Breast": "928286522",
              "Eye": "44026578",
              "Fallopian_Tube": "1739942440",
              "Heart": "1759721736",
              "Kidney": "949267305",
              "Knee": "1824489301",
              "Large_Intestine": "2043181688",
              "Liver": "1460762432",
              "Lung": "1824552484",
              "Lymph_Node": "1223566381",
              "Lymph_Vasculature": "1700987638",
              "Ovary": "1997082517",
              "Pancreas": "801179416",
              "Peripheral_Nervous_System": "714133140",
              "Placenta": "231591207",
              "Prostate": "1757780481",
              "Skin": "269383687",
              "Small_Intestine": "1762589435",
              "Spleen": "69626346",
              "Spinal_Cord": "1106564583",
              "Thymus": "863370556",
              "Ureter": "1106564583",
              "Urinary_Bladder": "1342577957",
              "Uterus": "1434605386"
             }


# In[15]:


version_dict = {"Blood": "v1.1",
                "Blood_Vascular": "v1.1",
                "Bone_Marrow": "v1.1",
                "Brain": "v1.1",
                "Breast": "v1.0",
                "Eye": "v1.0",
                "Fallopian_Tube": "v1.0",
                "Heart": "v1.1",
                "Kidney": "v1.1",
                "Knee": "v1.0",
                "Large_Intestine": "v1.1",
                "Liver": "v1.0",
                "Lung": "v1.1",
                "Lymph_Node": "v1.1",
                "Lymph_Vasculature": "v1.0",
                "Ovary": "v1.0",
                "Pancreas": "v1.0",
                "Peripheral_Nervous_System": "v1.0",
                "Placenta": "v1.0",
                "Prostate": "v1.0",
                "Skin": "v1.1",
                "Small_Intestine": "v1.0",
                "Spleen": "v1.1",
                "Spinal_Cord": "v1.0",
                "Thymus": "v1.1",
                "Ureter": "v1.0",
                "Urinary_Bladder": "v1.0",
                "Uterus": "v1.0"
               }


# ### Write ASCT+B URL to text file

# In[16]:


organ_list = [organ for organ, id in sheetID_dict.items()]
url_list = list()

read_url = open("url.txt", "w")
for organ in organ_list:  
    url = "https://docs.google.com/spreadsheets/d/{}/edit#gid={}".format(sheetID_dict[organ], gidID_dict[organ])
    url_list.append(url)
    read_url.write("Organ: {} - {}\n".format(organ, url))

#pprint.pprint(url_list)


# In[ ]:





# In[ ]:




