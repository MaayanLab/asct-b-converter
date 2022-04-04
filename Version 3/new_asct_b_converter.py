
# Harmonize ETL: HuBMAP
# Data Source: https://hubmapconsortium.github.io/ccf-asct-reporter/vis?selectedOrgans=All_Organs-v1.1,Blood-v1.1,Blood_Vasculature-v1.1,bone_marrow-v1.1,brain-v1.1,Eye-v1.0,Fallopian_Tube-v1.1,heart-v1.1,kidney-v1.1,Knee-v1.0,large_intestine-v1.1,Liver-v1.0,lung-v1.1,lymph_nodes-v1.1,Lymph_Vasculature-v1.0,Ovary-v1.0,Pancreas-v1.0,Peripheral_Nervous_System-v1.0,Prostate-v1.0,skin-v1.1,Small_Intestine-v1.0,spleen-v1.1,thymus-v1.1,Ureter-v1.0,Urinary_Bladder-v1.0,Uterus-v1.0&playground=false

import re
import sys
import os
import json
import requests
from bs4 import BeautifulSoup
#import pprint 

# # Initialization
# Input files should be named organ_name.tsv

# Enter input file name
file_name = sys.argv[1]
print("Working on file_name {}".format(file_name))
read_file = open(file_name, "r")

# Create a dictionary where the keys are columns (CT/BGene) and values are rows
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

# Output of column/row dictionary
excel_dict, key_list, version = tsv_to_columns_row_dictionary(read_file)

def get_jdata(url):
    myResponse = requests.get(url)
    if (myResponse.ok):
        j_data = json.loads(myResponse.content.decode())
        cell_type = j_data['_embedded']['terms'][0]['label']

        # Case where the cell ontology term is replaced 
        if "CL" in cell_type:
            cell_id = j_data['_embedded']['terms'][0]['term_replaced_by']
            url = "http://www.ebi.ac.uk/ols/api/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F{}".format(cell_id)
            cell_type = get_jdata(url)

    return cell_type


# Create a list of cell types. Cell types are listed at the highest resolution base on ASCT+B Table
# Key_list is a list of column name. Ex. CT/2... 
# Add cell type to dictionary as key
def add_cell_type_to_list(excel_dict, key_list):
    cell_type_list = list()
    as_name = excel_dict["AS/1"][0]
    for index in range(0, len(excel_dict["AS/1"])):
        cell_type_list.append(list())
    
    for index in range(0, len(excel_dict["AS/1"])):
        for label in key_list:
            if "CT" in label and "LABEL" not in label and "ID" not in label and "NOTES" not in label and "Reference" not in label and "DOI" not in label:
                if excel_dict[label][index]:
                    cell_type_list[index] = (excel_dict[label][index])
            if "ID" in label and "CT" in label and excel_dict[label][index]:
                try:
                    cell_id = excel_dict[label][index].replace(":", "_") 
                    url = "http://www.ebi.ac.uk/ols/api/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F{}".format(cell_id)
                    cell_type = get_jdata(url)
                    pprint.pprint(cell_type)
                    cell_type_list[index] = cell_type
                except:
                    print("Fetch cell type error")
                    continue
    
    # Dont add unique identifier
    final_cell_type_list = list()
    for cell_type in cell_type_list:
        final_cell_type_list.append(cell_type)
    
    return(as_name, final_cell_type_list)

# Cell Types List
as_name, cell_type_list = add_cell_type_to_list(excel_dict, key_list)
#pprint.pprint(cell_type_list)


def get_url(url):
    myResponse = requests.get(url)
    if (myResponse.ok):
        try:
            xml_data = myResponse.content
            soup = BeautifulSoup(xml_data, "xml")
            result_tag = soup.find("result")
            gene_tag = result_tag.find_all('str')
            gene_symbol = gene_tag[1].get_text()
        except:
            gene_symbol = ""

    else:
        gene_symbol = ""
    return gene_symbol

# Create a dictionary where the cell type is the key and the value is a gene marker
def add_gene_marker(cell_type_list, excel_dict):
    cell_type_gene_marker_dict = dict()
    
    for count, cell_type in enumerate(cell_type_list):
        if cell_type:
            if cell_type in cell_type_gene_marker_dict:   
                cell_type_gene_marker_dict[cell_type].update(excel_dict["All Gene Biomarkers"][count].split(", "))            
            else:
                cell_type_gene_marker_dict[cell_type] = set(excel_dict["All Gene Biomarkers"][count].split(", "))
    
    return(cell_type_gene_marker_dict)


def add_gene_protein_marker(cell_type_list, excel_dict, key_list):
    cell_type_gene_marker_dict = dict()
    
    for count, cell_type in enumerate(cell_type_list):
        if cell_type and excel_dict["Gene Markers"][count]:
            if cell_type in cell_type_gene_marker_dict:   
                cell_type_gene_marker_dict[cell_type].update(excel_dict["Gene Markers"][count].split(", "))            
            else:
                cell_type_gene_marker_dict[cell_type] = set(excel_dict["Gene Markers"][count].split(", "))

    # Add protein markers
    gene_list = list()
    for index in range(0, len(excel_dict["AS/1"])):
        gene_list.append(list())
        for label in key_list:
            if "BProtein" in label and "ID" in label and "LABEL" not in label and "Reference" not in label and "DOI" not in label:
                cell_value = excel_dict[label][index]
                cell_value_list = re.findall('HGNC:\d{4,5}', cell_value)
                cell_value_list = [i for i in cell_value_list if i]

                for cell_value in cell_value_list:
                    if cell_value:
                        url = "https://rest.genenames.org/search/hgnc_id/{}".format(cell_value.split(":")[1])
                        marker_symbol = get_url(url)
                        if marker_symbol:                    
                            gene_list[index].append(marker_symbol)

    for count, cell_type in enumerate(cell_type_list):
        if cell_type and gene_list[count]:
            if cell_type not in cell_type_gene_marker_dict.keys():
                cell_type_gene_marker_dict[cell_type] = set(gene_list[count])      
            else:
                cell_type_gene_marker_dict[cell_type].update(gene_list[count])
    
    return(cell_type_gene_marker_dict)

def add_gene_marker_no_all(cell_type_list, excel_dict, key_list):
    gene_list = list()
    for index in range(0, len(excel_dict["AS/1"])):
        gene_list.append(list())
        for label in key_list:
            cell_value = excel_dict[label][index]
            cell_value_list = re.findall('HGNC:\d{4,5}', cell_value)
            cell_value_list = [i for i in cell_value_list if i]

            for cell_value in cell_value_list:
                if cell_value:
                    url = "https://rest.genenames.org/search/hgnc_id/{}".format(cell_value.split(":")[1])
                    marker_symbol = get_url(url)
                    if marker_symbol: 
                        gene_list[index].append(marker_symbol)

    cell_type_gene_marker_dict = dict()
    for count, cell_type in enumerate(cell_type_list):
        if cell_type:
            if cell_type not in cell_type_gene_marker_dict.keys():
                cell_type_gene_marker_dict[cell_type] = set(gene_list[count])      
            else:
                cell_type_gene_marker_dict[cell_type].update(gene_list[count])
    
    return(cell_type_gene_marker_dict)


# Output of cell type/gene marker dictionary
#final_dict = add_gene_marker(cell_type_list, excel_dict, key_list)
# Output of cell type/gene marker dictionary
if "All Gene Biomarkers" in excel_dict:
    final_dict = add_gene_marker(cell_type_list, excel_dict)
elif "Gene Markers" in excel_dict:
    final_dict = add_gene_protein_marker(cell_type_list, excel_dict, key_list)
else:
    final_dict = add_gene_marker_no_all(cell_type_list, excel_dict, key_list)
#pprint.pprint(final_dict)

# Check if output directory exist
output_path = '{}/outdir_3'.format(os.getcwd())
if os.path.isdir(output_path):
    # Change the current working directory to output directory
    os.chdir(output_path)
else:
    # Create output directory if it does not exist
    os.mkdir(output_path)
    os.chdir(output_path)

output_filename = "{}_GMT_{}.tsv".format(file_name.split("/")[-1][:-4], version)
fh = open(output_filename , "w")
for key, value in final_dict.items():
    key = key.strip('\t') + ':' + as_name.title()
    if "[" in key or not value:
        continue
    elif value and type(value) == list:
        fh.write("{}\t{}\n".format(key.strip('\t'), ';'.join(value)))
    elif value:
        fh.write("{}\t{}\n".format(key, ';'.join(value)))
    else:
        fh.write("{}\n".format(key))



