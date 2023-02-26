scalar jDissPower(0);
scalar jFlowUni(0);


for (iLoop=0; iLoop<nObjPatch; iLoop++)
{
	if (objFunction==1) {
		jDissPower = jDissPower - sum(phi.boundaryField()[objPatchList[iLoop]]*(p.boundaryField()[objPatchList[iLoop]] + 0.5* magSqr(U.boundaryField()[objPatchList[iLoop]])));
		}
	else if (objFunction==2) {
		jFlowUni = jFlowUni - sum(0.5* magSqr(U.boundaryField()[objPatchList[iLoop]]));
		}
}

	if (objFunction==1) {
		Info<<"Objective Function (Power Dissipated) J: " <<jDissPower<<endl;
		}
	else if (objFunction==2) {
		Info<<"Objective Function (Flow uni @ outlet) J: " <<jFlowUni<<endl;
		}


