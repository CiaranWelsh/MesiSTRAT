import os, glob
import pandas, numpy
from collections import OrderedDict
import tellurium as te
import site
site.addsitedir(r'D:\pycotools3')
from pycotools3 import model, tasks, viz

# WORKING_DIRECTORY = os.path.join(os.path.dirname(__file__))
WORKING_DIRECTORY = r'D:\MesiSTRAT\BreastCancer'
MODELS_DIRECTORY = os.path.join(WORKING_DIRECTORY, 'models')
DATA_DIRECTORY = os.path.join(WORKING_DIRECTORY, 'data')
SS_DATA_FILE = fname = os.path.join(DATA_DIRECTORY, 'ss_data.csv')


COPASI_FILE = os.path.join(MODELS_DIRECTORY, 'simple_akt_model.cps')
COPASI_FORMATTED_DATA = fname = os.path.join(DATA_DIRECTORY, 'copasi_formatted_data.csv')

model_string = """
model SimpleAktModel()
    compartment Cell = 1;
    var IRS1    in Cell;
    var IRS1pS636_639   in Cell;
    var Akt     in Cell;
    var AktpT308    in Cell;
    var TSC2    in Cell;
    var TSC2pT1462   in Cell;
    var PRAS40  in Cell;
    var PRAS40pS183 in Cell;
    
    var FourEBP1    in Cell;
    var FourE_BP1pT37_46   in Cell;    
    const Insulin in Cell;
    
    // global variables
    Insulin = 0.005;
    // IRS1_obs                    := IRS1      + IRS1pS636_639   
    IRS1pS636_639_obs           := IRS1pS636_639          
    // Akt_obs                     := Akt       + AktpT308     
    AktpT308_obs                := AktpT308           
    // TSC2_obs                    := TSC2      + TSC2pT1462               
    TSC2pT1462_obs              := TSC2pT1462              
    // PRAS40_obs                  := PRAS40    + PRAS40pS183
    PRAS40pS183_obs             := PRAS40pS183        
    // FourEBP1_obs                 := FourEBP1   + FourE_BP1pT37_46
    // FourE_BP1pT37_46_obs        := FourE_BP1pT37_46          
    // S6K_obs                     := S6K       + S6KpT389
    S6KpT389_obs                := S6KpT389    
    
    //initial conditions
    // may need to do the 'total math'
    IRS1                = 1.700789 - IRS1pS636_639;
    IRS1pS636_639       = 0.861333;
    Akt                 = 1.241997 - AktpT308;
    AktpT308            = 0.486243;
    TSC2                = 1.136033 - TSC2pT1462;
    TSC2pT1462          = 0.644957;
    PRAS40              = 0.981968 - PRAS40pS183;
    PRAS40pS183         = 0.516932;
    // FourEBP1            = 0.458272 - FourE_BP1pT37_46;
    // FourE_BP1pT37_46    = 0.488169;
    S6K                 = 1.330735 - S6KpT389;
    S6KpT389            = 0.395656;
    
    // kinetic parameters
    _kIRS1Phos           = 0.1;     
    _kIRS1Dephos         = 0.1;         
    _kAktPhos            = 0.1;     
    _kAktDephos          = 0.1;     
    _kTSC2Phos           = 0.1;     
    _kTSC2Dephos         = 0.1;         
    _kPras40PhosByAkt    = 0.1;             
    _kPras40PhosByTSC    = 0.1;             
    _kPras40Dephos       = 0.1;         
    _kFourEBP1Phos           = 0.1;     
    _kFourEBP1Dephos         = 0.1;         
    _kS6KPhos            = 0.1;     
    _kS6KDephos          = 0.1;     
    
    // reactions
    R1 : IRS1 => IRS1pS636_639 ;        Cell*   _kIRS1Phos*Insulin*IRS1;
    R2 : IRS1pS636_639 => IRS1 ;        Cell*   _kIRS1Dephos*S6KpT389*IRS1pS636_639;
    R3 : Akt => AktpT308 ;          Cell*   _kAktPhos*IRS1pS636_639*Akt;
    R4 : AktpT308 => Akt ;          Cell*   _kAktDephos*AktpT308;
    R5 : TSC2 => TSC2pT1462 ;        Cell*   _kTSC2Phos*AktpT308*TSC2;
    R6 : TSC2pT1462 => TSC2 ;        Cell*   _kTSC2Dephos*TSC2pT1462;
    R7 : PRAS40 => PRAS40pS183 ;    Cell*   _kPras40PhosByAkt*AktpT308*PRAS40;
    R8 : PRAS40 => PRAS40pS183 ;    Cell*   _kPras40PhosByTSC*TSC2pT1462*PRAS40;
    R9 : PRAS40pS183 => PRAS40 ;    Cell*   _kPras40Dephos*TSC2pT1462;
    // R10: FourEBP1 => FourE_BP1pT37_46 ;        Cell*   _kFourEBP1Phos*TSC2*FourEBP1;
    // R11: FourE_BP1pT37_46 => FourEBP1 ;        Cell*   _kFourEBP1Dephos*FourE_BP1pT37_46;
    R12: S6K => S6KpT389 ;          Cell*   _kS6KPhos*TSC2*S6K;
    R13: S6KpT389 => S6K ;          Cell*   _kS6KDephos*S6KpT389;

end


"""

# model = te.loada(model_string)

mod = model.loada(model_string, copasi_file=COPASI_FILE)

with tasks.ParameterEstimation.Context(mod, [COPASI_FORMATTED_DATA, SS_DATA_FILE ],
                                       parameters='g') as context:
    context.set('separator', '\t')
    context.set('run_mode', False)
    context.set('randomize_start_values', True)
    context.set('method', 'genetic_algorithm')
    context.set('population_size', 100)
    context.set('lower_bound', 1e-1)
    context.set('upper_bound', 1e1)
    context.set('weight_method', 'standard_deviation')
    context.set('prefix', '_')
    config = context.get_config()

pe = tasks.ParameterEstimation(config)
mod = pe.models['simple_akt_model'].model



print(mod.open())

'''

Build a way to 'duplicate for each experiment' with pycotools. 
Set the initial concentrations to actual values during the 
parameter estimation by having lower and upper bounds equal to 
the correct value. 

'''






