{
      // scalarField sigmaMin (pow(mesh.V(),-1.0/3.0));
      //scalarField  toto ("toto",dimensionSet(0,-1,0,0,0,0,0));
      
   /*  Info << "ELSADX =" << ELSADx << endl;
     Info << "ELSADX.dimensions =" << ELSADx.dimensions() << endl;
     Info << "min sigmaMin =" << max(sigmaMin)<< endl;
     Info << "mesh.V() dimensions =" << mesh.V().dimensions() << endl;
     Info << "sigmaPrime dimensions =" << sigmaPrime.dimensions() << endl;

     Info << "ELSASigmaMin. dimensions =" << ELSASigmaMin.dimensions() << endl;
     */

    // alpha1.internalField()=max(alpha1.internalField(),0);
    // alpha1.internalField()=min(alpha1.internalField(),1);

     ELSARho=rho;

     Info << "compute ELSA SigmaMin ..." << endl;

     ELSASigmaMin.internalField()=CstELSASigmaMin*sqrt(alpha1*(1.-alpha1))*pow(mesh.V(),-1.0/3.0);


     Info << "compute ELSA SigmaEquil ..." << endl;

     if ( webEquilMod.value() < 0.5)
     {
       // WeberEquil Lebas 2009 IJMF
       ELSASigmaEquil.internalField()=rho*turbulence->k()/ELSAWeberEquil/
                                    ELSASurfaceTensionCoef*alpha1*(1.-alpha1);
       Info << "mod Lebas" << endl;
     }
     else
     {
       // Modified WeberEquil Benjamin Duret et al IJMF 2013
       ELSASigmaEquil.internalField()=2.*(rho1+rho2)*turbulence->k()/ELSAWeberEquil/
                                    ELSASurfaceTensionCoef*alpha1*(1.-alpha1);
       Info << "mod Duret "  << endl;
     }



     ELSASigmaEquil.max(0.0);

//    Not completely satisfactory, but it is to go over the undertermination of 
//    SigmaMin that goes to infinity as DX tends to zerO  pb in average and in 2D case ...

       ELSASigmaMin= min(ELSASigmaMin,ELSASigmaEquil*0.001);

//     Info << "compute ELSA Sigma ..." << endl;
 
 //    Info << ELSASigmaMin.boundaryField() << endl; 
 //    Info << ELSASigmaPrime.boundaryField() << endl; 

     ELSASigma.internalField()=ELSASigmaPrime.internalField()+ELSASigmaMin.internalField();

//     Info << ELSASigma.boundaryField() << endl; 

     ELSASigma.internalField()=max(ELSASigma.internalField(),0.0);


//   Info << ELSASigma.boundaryField() << endl; 


//    Info << "compute ELSA Ratio ..." << endl;

 //    ELSASigmaRatio=(turbulence->epsilon()/(turbulence->k()+kTiny))*(ELSASigma/(ELSASigmaEquil+sigmaTiny));
     ELSASigmaRatio.internalField()=(max((ELSASigma/(ELSASigmaEquil+sigmaTiny)),1.));
     ELSASigmaRatio.boundaryField()=ELSASigmaPrime.boundaryField();

//     Info << ELSASigmaRatio.boundaryField() << endl; 

//    Info << "... done" << endl;

     fvScalarMatrix ELSASigmaPrimeEqn
    (
        fvm::ddt(ELSASigmaPrime)
      + fvm::div(phi,ELSASigmaPrime)
      - fvm::laplacian
        (
            Dab + alphatab*turbulence->nut(), ELSASigmaPrime,
            "laplacian(Dab,ELSASigmaPrime)"
        )       
      == (CTurbELSASigmaPrime* (turbulence->epsilon()/(turbulence->k()+kTiny))  *(ELSASigmaRatio)  )*
             ( ELSASigmaEquil - (ELSASigmaPrime+ELSASigmaMin) )
//              (ELSASigmaPrime+ELSASigmaMin)*( 1.0 - (ELSASigma)/ (ELSASigmaEquil+sigmaTiny))
            
    );

/*    fvScalarMatrix ELSASigmaPrimeEqn
    (
        fvm::ddt(ELSASigmaPrime)
      + fvm::div(phi,ELSASigmaPrime)
      - fvm::laplacian
        (
            Dab + alphatab*turbulence->nut(), ELSASigmaPrime,
            "laplacian(Dab,ELSASigmaPrime)"
        )       
      == 0.2* ( ELSASigmaEquil - (ELSASigmaPrime+ELSASigmaMin) )/runTime.deltaT()
    ); */

    Info << "Solve ELSASigmaPrime ..." << endl;
    ELSASigmaPrimeEqn.relax();
    ELSASigmaPrimeEqn.solve();

    ELSASigmaPrime.max(0.0);

    ELSASigma=ELSASigmaPrime+ELSASigmaMin;
    Info << "...done" <<endl;
// D32=6*alpha/Sigma = 6*alpha/(SigmaMin+SigmaPrime)
// Sigmamin => Surface of a plane surface in mesh_vol_Cell
// D32 : a unique size of droplet or buble alpha=n*pi*D^3/6 and Sigma=n * pi*D^2
// if sigmaPrime is small it means that the previous definition i ....s false since the surface is plan 
// at mesh cell size order thus its curvature is more than delta x
// Dsigma = f(alpha)/(g(sigmamin)+h(sigmaprime)) = dx if sigmaprime --> 0 and equal 
// to 6*alpha/Sigma for sigmaprime large enough
   
// Dsigma(sigmaprime)=Dx-sigmaprime*f(alpha)/(g(sigmamin)+h(sigmaprime))^2*h'(sigmaprime)

    ELSAD32=6.0*alpha1*(1.0-alpha1)/(ELSASigmaMin+ELSASigmaPrime+ELSATiny);
}
