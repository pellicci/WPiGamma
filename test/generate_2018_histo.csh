#mkdir histos
#mkdir histos/latest_production_3referee

rm histos/latest_production/*2018*.root

#python generate_histos.py BDT 2 crab_projects/samples_MC_2018_3referee/latest_production/MC/backgrounds/WPiGammaAnalysis_DY10to50_2018.root histos/latest_production_3referee/WPiGammaHistos_DY10to50_2018.root
#python generate_histos.py BDT 2 crab_projects/samples_MC_2018_3referee/latest_production/MC/backgrounds/WPiGammaAnalysis_DY50_2018.root histos/latest_production_3referee/WPiGammaHistos_DY50_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_DY10to50_2018.root histos/latest_production/WPiGammaHistos_DY10to50_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_DY50_2018.root histos/latest_production/WPiGammaHistos_DY50_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_GammaJets20to40_2018.root histos/latest_production/WPiGammaHistos_GammaJets20to40_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_GammaJets20toInf_2018.root histos/latest_production/WPiGammaHistos_GammaJets20toInf_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_GammaJets40toInf_2018.root histos/latest_production/WPiGammaHistos_GammaJets40toInf_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_QCDDoubleEMEnriched30to40_2018.root histos/latest_production/WPiGammaHistos_QCDDoubleEMEnriched30to40_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_QCDDoubleEMEnriched30toInf_2018.root histos/latest_production/WPiGammaHistos_QCDDoubleEMEnriched30toInf_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_QCDDoubleEMEnriched40toInf_2018.root histos/latest_production/WPiGammaHistos_QCDDoubleEMEnriched40toInf_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_QCDHT100to200_2018.root histos/latest_production/WPiGammaHistos_QCDHT100to200_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_QCDHT200to300_2018.root histos/latest_production/WPiGammaHistos_QCDHT200to300_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_QCDHT300to500_2018.root histos/latest_production/WPiGammaHistos_QCDHT300to500_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_QCDHT500to700_2018.root histos/latest_production/WPiGammaHistos_QCDHT500to700_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_QCDHT700to1000_2018.root histos/latest_production/WPiGammaHistos_QCDHT700to1000_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_QCDHT1000to1500_2018.root histos/latest_production/WPiGammaHistos_QCDHT1000to1500_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_QCDHT1500to2000_2018.root histos/latest_production/WPiGammaHistos_QCDHT1500to2000_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_QCDHT2000toInf_2018.root histos/latest_production/WPiGammaHistos_QCDHT2000toInf_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_SingleAntiToptW_2018.root histos/latest_production/WPiGammaHistos_SingleAntiToptW_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_SingleToptW_2018.root histos/latest_production/WPiGammaHistos_SingleToptW_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_TTGJets_2018.root histos/latest_production/WPiGammaHistos_TTGJets_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_WGToLNuG01J_2018.root histos/latest_production/WPiGammaHistos_WGToLNuG_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_WJetsToLNu0J_2018.root histos/latest_production/WPiGammaHistos_WJetsToLNu0J_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_WJetsToLNu1J_2018.root histos/latest_production/WPiGammaHistos_WJetsToLNu1J_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_WJetsToLNu2J_2018.root histos/latest_production/WPiGammaHistos_WJetsToLNu2J_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_WW_2018.root histos/latest_production/WPiGammaHistos_WW_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_WZ_2018.root histos/latest_production/WPiGammaHistos_WZ_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_ZGTo2LG_2018.root histos/latest_production/WPiGammaHistos_ZGTo2LG_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_ttbarWQQ_2018.root histos/latest_production/WPiGammaHistos_ttbarWQQ_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_ttbarWlnu_2018.root histos/latest_production/WPiGammaHistos_ttbarWlnu_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_ttbarZQQ_2018.root histos/latest_production/WPiGammaHistos_ttbarZQQ_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_ttbarZlnu_2018.root histos/latest_production/WPiGammaHistos_ttbarZlnu_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_ttbarlnu_2018.root histos/latest_production/WPiGammaHistos_ttbarlnu_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_ttbarToHadronic_2018.root histos/latest_production/WPiGammaHistos_ttbarToHadronic_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/backgrounds/WPiGammaAnalysis_ttbarToSemiLeptonic_2018.root histos/latest_production/WPiGammaHistos_ttbarToSemiLeptonic_2018.root

#python generate_histos.py BDT 2 crab_projects/samples_MC_2018_3referee/latest_production/MC/signals/WPiGammaAnalysis_Signal_2018.root histos/latest_production_3referee/WPiGammaHistos_Signal_2018.root
python generate_histos.py BDT 2 rootfiles/latest_production/MC/signals/WPiGammaAnalysis_Signal_2018.root histos/latest_production/WPiGammaHistos_Signal_2018.root

python generate_histos.py BDT 2 rootfiles/latest_production/dataprocess/WPiGammaAnalysis_Data_2018.root histos/latest_production/WPiGammaHistos_Data_2018.root

#Merge samples
hadd histos/latest_production/WPiGammaHistos_QCD_2018.root histos/latest_production/WPiGammaHistos_QCDHT*_2018.root histos/latest_production/WPiGammaHistos_QCDDouble*_2018.root
#hadd histos/latest_production/WPiGammaHistos_QCDEM_2018.root histos/latest_production/WPiGammaHistos_QCDDouble*_2018.root
rm histos/latest_production/WPiGammaHistos_QCDHT*_2018.root
rm histos/latest_production/WPiGammaHistos_QCDDouble*_2018.root

hadd histos/latest_production/WPiGammaHistos_GammaJets_2018.root histos/latest_production/WPiGammaHistos_GammaJets*to*_2018.root
rm histos/latest_production/WPiGammaHistos_GammaJets*to*_2018.root

hadd histos/latest_production/WPiGammaHistos_STtW_2018.root histos/latest_production/WPiGammaHistos_Single*ToptW_2018.root
rm histos/latest_production/WPiGammaHistos_Single*ToptW_2018.root

hadd histos/latest_production/WPiGammaHistos_DY_2018.root histos/latest_production/WPiGammaHistos_DY*50_2018.root
rm histos/latest_production/WPiGammaHistos_DY*50_2018.root

hadd histos/latest_production/WPiGammaHistos_ttbar_2018.root histos/latest_production/WPiGammaHistos_ttbarTo*nic_2018.root histos/latest_production/WPiGammaHistos_ttbarlnu_2018.root
rm histos/latest_production/WPiGammaHistos_ttbarTo*nic_2018.root
rm histos/latest_production/WPiGammaHistos_ttbarlnu_2018.root

hadd histos/latest_production/WPiGammaHistos_WJetsToLNu_2018.root histos/latest_production/WPiGammaHistos_WJetsToLNu*J_2018.root
rm histos/latest_production/WPiGammaHistos_WJetsToLNu*J_2018.root
