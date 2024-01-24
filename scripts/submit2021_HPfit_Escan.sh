### ENERGY SCAN ON C3 ###
# with LP template @ 100 GeV
## 100 GeV
#python scripts/submitBatch.py -r 14918 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_LP-HP.cfg -v fitLP -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 
### 150 GeV
#python scripts/submitBatch.py -r 14934 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_LP-HP.cfg -v fitLP -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 
#python scripts/submitBatch.py -r 14943 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_LP-HP.cfg -v fitLP -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 
### 200 GeV
#python scripts/submitBatch.py -r 14951 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_LP-HP.cfg -v fitLP -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 
### 250 GeV
#python scripts/submitBatch.py -r 14820 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_LP-HP.cfg -v fitLP -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 
#python scripts/submitBatch.py -r 14821 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_LP-HP.cfg -v fitLP -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 

# with HP template x VFEs @ 100 GeV
##  50 GeV - merge run taken at the edge of the C3 crystal
python scripts/submitBatch.py -r 14907 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow # U 
python scripts/submitBatch.py -r 14914 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow # D 
python scripts/submitBatch.py -r 14931 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow # L
python scripts/submitBatch.py -r 14932 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow # R
python scripts/submitBatch.py -r 14937 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow # R
python scripts/submitBatch.py -r 14938 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow # L
## 100 GeV
python scripts/submitBatch.py -r 14918 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 
## 150 GeV
python scripts/submitBatch.py -r 14934 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 
python scripts/submitBatch.py -r 14943 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 
python scripts/submitBatch.py -r 14940 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow #U 
python scripts/submitBatch.py -r 14942 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow #R 
python scripts/submitBatch.py -r 14944 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow #L 
python scripts/submitBatch.py -r 14949 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow #D 
## 200 GeV
python scripts/submitBatch.py -r 14951 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 
## 250 GeV
python scripts/submitBatch.py -r 14820 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 
python scripts/submitBatch.py -r 14821 -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v fitVFEs_fixes -s /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2021/HighPurity/ --notar  --spills-per-job 1 -q tomorrow 


