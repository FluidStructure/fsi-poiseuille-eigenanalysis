//==============|
//  NAME        : kernellibclass.cpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 11.10.2008
//  DESCRIPTION : Holds the kernel functions for many different particle interactions.
//  NOTES       : Suitable for use with the fmm or as a standalone function.
//  TODO        : .
//==============|
#include "kernellib.hpp"

//  ------------:
//  NAME        : lambVortex
//  PURPOSE     : This is a solver for the velocity induced by a discrete lamb vortex on a point in 2d space.
//  IMPORTS     : .
//
//  PRE-CONDS   : .
//  POST-CONDS  : .
//  ------------:
ComplexDouble KernelLib::lambVortex(ComplexDouble particlePos, double particleStr, double coreSizeSqrd, ComplexDouble evalPos, double minPrecision)
{
    ComplexDouble r;
    double rSquared, u,v, vel, rCuttoffSqrd, coreCorr, ZERO;

    r = evalPos-particlePos;
    rSquared = norm(r);
    vel = 0.0;
    ZERO = 1e-16; // This could be a problem? could make this relative to core size to an acceptable precision?
    if (rSquared > ZERO)
    {
        // presuming particlestr is already divided by 2*pi
        vel = particleStr / rSquared;

        if (minPrecision > 0.0)
        {
            rCuttoffSqrd = coreSizeSqrd * -log(minPrecision); // point vortex assumption is true to tolerance only outside of this zone
        }
        else
        {
            // want exact precision so will use minPrecision <=0.0 condition to complete the calculation
            rCuttoffSqrd = 0.0;
        }

        if ((rSquared < rCuttoffSqrd) || (minPrecision <= 0.0))
        {
            //	Is a lamb vortex where the evaluation point is located inside
            //	the core
            coreCorr = 1.0 - exp(-1.0 * rSquared / (coreSizeSqrd));
            vel *= coreCorr;
        }
    }
    u=vel*imag(r);
    v=-vel*real(r);

    return ComplexDouble(u,v);    //    return u and v velocity values through outVel
}

//  ------------:
//  NAME        : lambVortexVort
//  PURPOSE     : This is a solver for the vorticity induced by a discrete lamb vortex on a point in 2d space.
//  IMPORTS     : .
//
//  PRE-CONDS   : .
//  POST-CONDS  : .
//  ------------:
ComplexDouble KernelLib::lambVortexVort(ComplexDouble particlePos, double particleStr, double coreSizeSqrd, ComplexDouble evalPos, double minPrecision)
{
    ComplexDouble r;
    double rSquared, vort, rCuttoffSqrd, ZERO;

    r = evalPos-particlePos;
    rSquared = norm(r);
    vort = 0.0;
    ZERO = 1e-16;    // This could be a problem? could make this relative to core size to an acceptable precision?
    if (coreSizeSqrd > ZERO)
    {
        rCuttoffSqrd = coreSizeSqrd * -log(minPrecision); // point vortex assumption is true to tolerance only outside of this zone
        if (rSquared < (rCuttoffSqrd))
        {
            //	Is a lamb vortex where the evaluation point is located inside
            //	the core according to prescribed precision
            vort = exp(-1.0 * rSquared / (coreSizeSqrd)) * particleStr / (coreSizeSqrd);
        }
    }

    return ComplexDouble(vort, 0.0);    //    return vorticity through a complex number
}
//  ------------:
//  NAME        : vortexElement
//  PURPOSE     : This is a solver for the velocity induced by a constant/linear strength vortex elements on a point in 2d space.
//  IMPORTS     : .
//
//  PRE-CONDS   : .
//  POST-CONDS  : .
//  NOTES       : Positive vorticity is clockwise
//  ------------:
ComplexDouble KernelLib::vortexElement(ComplexDouble posL, double strL, ComplexDouble posR, double strR, ComplexDouble evalPos, double panelTolerance, bool isInfiniteSheet)
{
    double vVel, uVel, ZERO, vP, uP, cosAlpha, sinAlpha;
    //for linear sheet, vel == vel from constant + so if vel from linear node1 != node2, gradient == node1-node2/length
    uVel = 0.0;
    vVel = 0.0;
    vP = 0.0;
    uP = 0.0;
    cosAlpha = 0.0;
    sinAlpha = 0.0;

    ComplexDouble panLR = posR - posL;
    double panLen = abs(panLR);
    if (!isInfiniteSheet)
    {
        ZERO = panLen * abs(panelTolerance);	//	1e-3 for tolerance good default
    }
    else
    {
        ZERO = 1e-16;       //not sure what other tolerance would be suitable for seminfinite!
    }

    if ((abs(posL - evalPos) > ZERO) && (abs(posR - evalPos) > ZERO))
    {
        //	Evaluation point is NOT on the end of the panel.
        ComplexDouble panCen = (posL + posR)/2.0;

        if (panLen > ZERO) // else it is a point vortex so callup lvtvel and display a warning msg
        {
            double strPerM, globalXRel, globalYrel;
            double xRel, yRel, xlRel, xrRel;
            double rRightSqrd, rLeftSqrd, log_Rr_Rl, thetaL, thetaR;
            ComplexDouble globalRel;

            strPerM = strL;

            // alpha in terms of katz and plotkins terminoligy
            cosAlpha =  panLR.real() / panLen;
            sinAlpha = -panLR.imag() / panLen;

            globalRel = evalPos - panCen;
            globalXRel = globalRel.real();
            globalYrel = globalRel.imag();

            //	convert to local panel coordinates relative to the panel!
            // Calculate the relative x and y distances from element "j"s center
            //   page 318 in katz and plotkins transformation matrices.
            xRel = globalXRel* cosAlpha - globalYrel* sinAlpha;
            yRel = globalXRel* sinAlpha + globalYrel* cosAlpha;
            if (abs(yRel) < ZERO)
            {
                // yRel is closer than the accepatble tolerance that we have set, hence reset the polarity and continue
                // cause if gmres grief when using python due to differing tolerances when talking about on center of the panel
                yRel = 0.0;
            }
            if (abs(xRel) < ZERO)
            {   // practically on centre of panel. ensures vertical component is zero
                xRel = 0.0;
            }       

            xlRel = xRel + (panLen/2.0);    // equivalent to k&p (x-x_1) terminology
            xrRel = xRel - (panLen/2.0);    // equivalent to k&p (x-x_2) terminology

            // These next terms relate directly to that you would find in the k&p equations
            thetaL = atan2(yRel, xlRel);
            thetaR = atan2(yRel, xrRel);
            rRightSqrd = pow(xrRel,2) + pow(yRel,2);
            rLeftSqrd = pow(xlRel,2) + pow(yRel,2);
            log_Rr_Rl = log(rRightSqrd/rLeftSqrd);       // this is repeated 3 times so may aswell save some cpu time! Never need to worry about 0/0 as thats ~condition of this current block

            // strength already been precorrected by /2PI
            uP = strPerM  * (thetaR - thetaL);    //   in k&p terms, see k&p eq-10.39 page.273 and/or eq10.18 p269
            vP = (strPerM / 2.0) * log_Rr_Rl;	          // see k&p eq-10.40 page.273 and/or eq10.17 p269

            if (isInfiniteSheet)
            {
                // these terms turn a panel velocity into a seminfinite sheet velocity, as a panel is just the hollow section of the infinite sheet
                uP = -uP + (strPerM * M_PI * sign(yRel));
                vP *= -1.0;
            }
            else if (abs(strR - strL) > ZERO)
            {
                // Constant panel component already done maybe we instead have a linear panel
                // do extra Linear panel calculations too
                // see k&p eq 10.72 page.279 and/or eq10.48 p269
                double strGradient = (strR - strL) / panLen;
                uP += strGradient * ( yRel * log_Rr_Rl / 2.0 + xlRel * (thetaR - thetaL));
                vP += strGradient * ( xlRel * log_Rr_Rl / 2.0 + yRel * (thetaL - thetaR) + panLen );
            }
        }
        else
        {
            // Panel has no length, so it must be a lamb vortex.
            cout << "Warning this panel has no length!, check your inputs are correct" << endl;
        }
    }

    // converting back to global coordinates
    uVel =   uP* cosAlpha + vP* sinAlpha;
    vVel = - uP* sinAlpha + vP* cosAlpha;

    return ComplexDouble(uVel,vVel);	//	return u and v velocity values
}