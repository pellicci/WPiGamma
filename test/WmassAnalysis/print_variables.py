import ROOT

#Get the model and the data
#fInput = ROOT.TFile("fitData.root")
fInput = ROOT.TFile("fitMC.root")
fInput.cd()

workspace = fInput.Get("workspace")

#workspace.Print()

argset =  workspace.allVars()

#v = workspace.var("a0_el")

_variable = argset.createIterator()
variable = _variable.Next()

while variable:

    print variable.GetName(), "            ", variable.getVal(), "          ", variable.isConstant()
    variable = _variable.Next()
