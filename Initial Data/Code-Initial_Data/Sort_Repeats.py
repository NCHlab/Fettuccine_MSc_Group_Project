#! /usr/bin/python
# coding: utf-8

"""This python script takes the initial tables exported by the Get_the_inital_tables.py script and creates new 
.csv files with sorted families, superfamilies, repeats"""
import pandas as pd
import re
    
def cat_by_families(Choice):
    
    HERV_reps = pd.read_csv("HERV_assign_to_families.csv")
        
    fams_names_group = HERV_reps.groupby(str(Choice))['Repeat Name'].unique()
    fams_IDs_group = HERV_reps.groupby(str(Choice))['ID'].unique()
    
    if Choice == "Family":
        fams_sfam_group = HERV_reps.groupby(str(Choice))['Superfamily'].unique()
        HERV_Fams_SFams = fams_sfam_group.to_frame()
  
    
    HERV_Fams_IDs = fams_IDs_group.to_frame()
    HERV_Fams_IDs["Count"] = HERV_Fams_IDs["ID"].apply(lambda x: x.size)

    HERV_Families = fams_names_group.to_frame()
    
    HERV_Families["ID"] = HERV_Fams_IDs.reindex(HERV_Families.index)["ID"]
    HERV_Families["Count"] = HERV_Fams_IDs.reindex(HERV_Families.index)["Count"]
    HERV_Families["Repeat Name"] = HERV_Families["Repeat Name"].apply(lambda x: ", ".join(x)) 
    HERV_Families["ID"] = HERV_Families["ID"].apply(lambda x: ", ".join(x)) 
    if Choice == "Family": 
        HERV_Families["Superfamily"] = HERV_Fams_SFams.reindex(HERV_Families.index)["Superfamily"]
        HERV_Families["Superfamily"] = HERV_Families["Superfamily"].apply(lambda x: ", ".join(x)) 
        
    HERV_Families.to_csv("HERV_"+str(Choice)+"_cat.csv")
cat_by_families("Family")
cat_by_families("Superfamily")


def RepeatsThemselves(Choice):
    if Choice == "HERV-U":
        HERV_reps = pd.read_csv("HERV_assign_to_families.csv")
        reps = HERV_reps[HERV_reps["Family"] == "Unassigned"]
        reps["Repeat Name"] = reps["Repeat Name"].apply(lambda x: x.replace("-int", ""))
        reps["Repeat Name"] = reps["Repeat Name"].apply(lambda x: x.replace("_$", "$"))
        reps["Repeat Name"] = reps["Repeat Name"].apply(lambda x: x.replace("-", ""))
        out = "HERV_Unassigned_reps_cat.csv"
    elif Choice == "HERV-A":
        HERV_reps = pd.read_csv("HERV_assign_to_families.csv")
        reps = HERV_reps[HERV_reps["Family"] != "Unassigned"]
        reps["Repeat Name"] = reps["Repeat Name"].apply(lambda x: x.replace("-int", ""))
        reps["Repeat Name"] = reps["Repeat Name"].apply(lambda x: x.replace("_$", "$"))
        reps["Repeat Name"] = reps["Repeat Name"].apply(lambda x: x.replace("-", ""))
        out = "HERV_Assigned_reps_cat.csv"
    elif Choice == "L1":
        reps = pd.read_csv("L1_names_family.csv")
        out = "L1_reps_cat.csv"
        
    else: print "Please type either Unassigned or Assigned"

    
    names_group = reps.groupby('Repeat Name')['ID'].unique()
    
    rep_cat = names_group.to_frame()
    rep_cat["Count"] = rep_cat["ID"].apply(lambda x: x.size)
    rep_cat["ID"] = rep_cat["ID"].apply(lambda x: ", ".join(x))
    rep_cat.to_csv(out)

print "What type of repeats? Type L1 or HERV:"
loh = raw_input(">>> ")
if loh == "L1":
    RepeatsThemselves("L1")
    print "L1 repeats sorted by themselves"
elif loh == "HERV":
    cat_by_families("Family")
    cat_by_families("Superfamily")
    print "HERVs sorted by Families and Superfamilies"
    print ""
    RepeatsThemselves("HERV-A")
    RepeatsThemselves("HERV-U")
    print "HERV repeats sorted by themselves"
    

    

