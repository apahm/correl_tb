
///////////////////////////////////////////////
//             Richard E. Bird               //
// Clear Sky Broadband Solar Radiation Model //
//                                           //
//            September 19, 2012             //
//                                           //
//   Filename: BIRD.C                        //
//                                           //
//   Afshin Michael Andreas                  //
//   Afshin.Andreas@NREL.gov (303)384-6383   //
//                                           //
//   Solar Resource and Forecasting Group    //
//   Solar Radiation Research Laboratory     //
//   National Renewable Energy Laboratory    //
//   15013 Denver W Pkwy, Golden, CO 80401   //
///////////////////////////////////////////////

///////////////////////////////////////////////
//   See the BIRD.H header file for usage    //
//                                           //
//  This code is based on the SERI (NREL)    //
//  technical report "A Simplified Clear     //
//  Sky model for Direct and Diffuse         //
//  Insolation on Horizontal Surfaces" by    //
//  R.E. Bird and R.L. Hulstrom              //
///////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
// Revised 11-MAY-2017 Andreas
//         Set DNI, GHI, or DHI to zero if the result from model is negative
//         Prevent posible divide by zero when calculating GHI
///////////////////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "bird.h"

void bird_horiz(double dni, double coszen, double ias, double albedo, double rs, double *ghi, double *dhi)
{
	double dir_horz = dni*coszen;

	if (albedo*rs != 1.0)
	{
		 *ghi = (dir_horz + ias) / (1.0 - albedo*rs);
		 if (*ghi < 0) *ghi = 0;
	} else *ghi = 0;

	*dhi = *ghi - dir_horz;
	if (*dhi < 0) *dhi = 0;
}

void bird_calculate(bird_data *b)
{
	double etrn, press, oz, wat, coszen, taa, rs, ias;
	double t_rayliegh, t_ozone, t_gases, t_water, t_aerosol;

	if ((b->zenith >= 0) && (b->zenith < 90) && (b->r > 0))
	{
		etrn     = 1367.0 / (b->r*b->r);
		coszen   = cos(b->zenith*3.14159265358979323846264338328/180.0);
		b->amass = 1.0 / (coszen + 0.50572 * pow((96.07995 - b->zenith),-1.6364));


		press  = b->pressure*b->amass/1013;
		oz     = b->ozone*b->amass;
		wat    = b->water*b->amass;

		t_rayliegh = exp(-0.0903*pow(press,0.84)*(1+press-pow(press,1.01)));
		t_ozone    = 1-0.1611*(oz)*pow(1+139.48*oz,-0.3034)-0.002715*(oz)/(1+0.044*(oz)+0.0003*oz*oz);
		t_gases    = exp(-0.0127*pow(press,0.26));
		t_water    = 1-2.4959*wat/(pow(1+79.034*wat,0.6828)+6.385*wat);
		t_aerosol  = exp(-(pow(b->taua,0.873))*(1+b->taua-pow(b->taua,0.7088))*pow(b->amass,0.9108));

		b->direct_normal = 0.9662*etrn*t_aerosol*t_water*t_gases*t_ozone*t_rayliegh;
		if (b->direct_normal < 0) b->direct_normal = 0;

		taa = 1-0.1*(1-b->amass+pow(b->amass,1.06))*(1-t_aerosol);
		rs  = 0.0685+(1-b->ba)*(1-t_aerosol/taa);
		ias = etrn*coszen*0.79*t_ozone*t_gases*t_water*taa*(0.5*(1-t_rayliegh)+b->ba*(1-(t_aerosol/taa)))/(1-b->amass+pow(b->amass,1.02));

		bird_horiz(b->direct_normal, coszen, ias, b->albedo, rs, &b->global_horiz, &b->diffuse_horiz);

		if (b->dni_mod >= 0)
		{
			b->direct_normal_mod = b->direct_normal * b->dni_mod;
			bird_horiz(b->direct_normal_mod, coszen, ias, b->albedo, rs, &b->global_horiz_mod, &b->diffuse_horiz_mod);
		}

	} else {
		b->amass = b->direct_normal = b->global_horiz = b->diffuse_horiz = 0;
		b->direct_normal_mod = b->global_horiz_mod = b->diffuse_horiz_mod = 0;
	}
}