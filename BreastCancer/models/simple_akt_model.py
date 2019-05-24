import os, glob
import pandas, numpy
from collections import OrderedDict
import tellurium as te
from pycotools3 import model, tasks, viz

model_string = """
model SimpleAktModel()
    compartment Cell = 1;
    var IRS1    in Cell;
    var pIRS1   in Cell;
    var Akt     in Cell;
    var pAkt    in Cell;
    var TSC1    in Cell;
    var pTSC1   in Cell;
    var PRAS40  in Cell;
    var pPRAS40 in Cell;
    
    var FourEBP    in Cell;
    var pFourEBP   in Cell;    
    const Insulin in Cell;
    
    // global variables
    Insulin = 1;
    
    //initial conditions
    IRS1    = 0.1;
    pIRS1   = 0.1;
    Akt     = 0.1;
    pAkt    = 0.1;
    TSC1    = 0.1;
    pTSC1   = 0.1;
    PRAS40  = 0.1;
    pPRAS40 = 0.1;
    FourEBP    = 0.1;
    pFourEBP   = 0.1;
    
    // kinetic parameters
    _kIRS1Phos           = 0.1;     
    _kIRS1Dephos         = 0.1;         
    _kAktPhos            = 0.1;     
    _kAktDephos          = 0.1;     
    _kTSC1Phos           = 0.1;     
    _kTSC1Dephos         = 0.1;         
    _kPras40PhosByAkt    = 0.1;             
    _kPras40PhosByTSC    = 0.1;             
    _kPras40Dephos       = 0.1;         
    _kFourEBPPhos           = 0.1;     
    _kFourEBPDephos         = 0.1;         
    _kS6KPhos            = 0.1;     
    _kS6KDephos          = 0.1;     
    
    // reactions
    R1 : IRS1 => pIRS1 ;        Cell*   _kIRS1Phos*Insulin*IRS1;
    R2 : pIRS1 => IRS1 ;        Cell*   _kIRS1Dephos*pS6K*pIRS1;
    R3 : Akt => pAkt ;          Cell*   _kAktPhos*pIRS1*Akt;
    R4 : pAkt => Akt ;          Cell*   _kAktDephos*pAkt;
    R5 : TSC1 => pTSC1 ;        Cell*   _kTSC1Phos*pAkt*TSC1;
    R6 : pTSC1 => TSC1 ;        Cell*   _kTSC1Dephos*pTSC1;
    R7 : PRAS40 => pPRAS40 ;    Cell*   _kPras40PhosByAkt*pAkt*PRAS40;
    R8 : PRAS40 => pPRAS40 ;    Cell*   _kPras40PhosByTSC*pTSC1*PRAS40;
    R9 : pPRAS40 => PRAS40 ;    Cell*   _kPras40Dephos*pTSC1;
    R10: FourEBP => pFourEBP ;        Cell*   _kFourEBPPhos*TSC1*FourEBP;
    R11: pFourEBP => FourEBP ;        Cell*   _kFourEBPDephos*pFourEBP;
    R12: S6K => pS6K ;          Cell*   _kS6KPhos*TSC1*S6K;
    R13: pS6K => S6K ;          Cell*   _kS6KDephos*pS6K;

end


"""

model = te.loada(model_string)

print(model)
