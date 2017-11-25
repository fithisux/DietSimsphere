/*
 * Copyright (C) 2015 Toby N. Carlson
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
 /*--- formatted by Jindent 2.1, (www.c-lab.de/~jindent) ---*/
package psu.ets.psubams.model;


/* Start of the new psubams architecture. */

 /*
 * The idea here is to recode all the fortran
 * as a single block of C, and then break it
 * down into functions afterwards. This means
 * that all the common variables can be common
 * for the time being, and the resulting functions
 * will (hopefully) be a bit more logical.
 */
 /* Standard header files */

 /*
 * #include <stdlib.h>
 * #include <stdio.h>
 * #include <h>
 */
 /* My header files */

 /*
 * #include "psubams.h"
 * #include "const.h"
 * #include "dectim.h"
 * #include "spline.h"
 * #include "rlinear.h"
 * #include "transm.h"
 * #include "velfns.h"
 * #include "stomfs.h"
 * #include "psgcal.h"
 * #include "stomc.h"
 * #include "daykm.h"
 * #include "fine.h"
 * #include "cond.h"
 * #include "misc.h"
 */
/**
 * Class declaration
 *
 *
 * @author
 * @version %I%, %G%
 */
public class Model implements Constants, SeriesNames {

    /**
     * Method declaration
     *
     *
     * @param iyr
     * @param imo
     * @param iday
     * @param tz
     * @param xlat
     * @param xlong
     * @param strtim
     * @param timend
     * @param outtt
     * @param slope
     * @param aspect
     * @param f
     * @param fsub
     * @param wmax
     * @param btemp
     * @param tp
     * @param dual_ti
     * @param ti_a
     * @param ti_b
     * @param albedo_gflag
     * @param albg
     * @param epsi
     * @param omega
     * @param zo
     * @param obst_hgt
     * @param cloud_flag
     * @param cld_fract
     * @param frveg
     * @param xlai
     * @param epsf
     * @param albedo_fflag
     * @param albf
     * @param stmtype
     * @param volrel
     * @param rmin
     * @param rcut
     * @param wilt
     * @param vegheight
     * @param width
     * @param steady
     * @param ci
     * @param co
     * @param coz_sfc
     * @param coz_air
     * @param rks
     * @param cosbyb
     * @param thmax
     * @param psis
     * @param mintemp
     * @param maxtemp
     * @param beta
     * @param b1
     * @param b2
     * @param psice
     * @param sc
     * @param zp
     * @param frhgt
     * @param frzp
     * @param rkocap
     * @param rccap
     * @param rzcap
     * @param volini
     * @param zstini
     * @param nobs_ptq
     * @param nobs_wind
     * @param station_height
     * @param ugs
     * @param vgs
     * @param ps
     * @param ts
     * @param dep
     * @param dd0
     * @param ff0
     * @param zh
     * @param output_data
     * @param sounding_data
     *
     * @throws Exception
     *
     * @see
     */
    public void _psubams(int iyr, int imo, int iday, double tz, double xlat,
            double xlong, double strtim, double timend,
            double outtt, double slope, double aspect, double f,
            double fsub, double wmax, double btemp, double tp,
            char dual_ti, double ti_a, double ti_b,
            char albedo_gflag, double albg, double epsi,
            double omega, double zo, double obst_hgt,
            char cloud_flag, double cld_fract, double frveg,
            double xlai, double epsf, char albedo_fflag,
            double albf, char stmtype, double volrel,
            double rmin, double rcut, double wilt,
            double vegheight, double width, char steady,
            double ci, double co, double coz_sfc,
            double coz_air, double rks, double cosbyb,
            double thmax, double psis, double mintemp,
            double maxtemp, double beta, double b1, double b2,
            double psice, double sc, double zp, double frhgt,
            double frzp, double rkocap, double rccap,
            double rzcap, double volini, double zstini,
            int nobs_ptq, int nobs_wind, double station_height,
            double ugs, double vgs, double[] ps, double[] ts,
            double[] dep, double[] dd0, double[] ff0,
            double[] zh,
            /* String[] output_names, */
            double[][] data_graphs, double[][][] sounding_data,
            int dayOfTheYear/*update by Yannis Konstas*/,
            double[][] output_data, boolean batch /*double[][][] sounding_data*/) throws Exception {

        /*	 System.out.println("frveg: " + frveg);
	System.out.println("rmin: " + rmin);
	System.out.println("b1: " + b1);
	System.out.println("iyr: " + iyr);
	System.out.println("imo: " + imo);
	System.out.println("iday: " + iday);
	System.out.println("tz: " + tz);
	System.out.println("xlat: " + xlat);
	System.out.println("xlong: " + xlong);
	System.out.println("strtim: " + strtim);
	System.out.println("timend: " + timend);
	System.out.println("outtt: " + outtt);
	System.out.println("slope: " + slope);
	System.out.println("aspect: " + aspect);
	System.out.println("f: " + f);
	System.out.println("fsub: " + fsub);
	System.out.println("wmax: " + wmax);
	System.out.println("btemp: " + btemp);
	System.out.println("tp: " + tp);
	System.out.println("dual_ti: " + dual_ti);
	System.out.println("ti_a: " + ti_a);
	System.out.println("ti_b: " + ti_b);
	System.out.println("albedo_gflag: " + albedo_gflag);
	System.out.println("albg: " + albg);
	System.out.println("epsi: " + epsi);
	System.out.println("omega: " + omega);
	System.out.println("zo: " + zo);
	System.out.println("obst_hgt: " + obst_hgt);
	System.out.println("cloud_flag: " + cloud_flag);
	System.out.println("cld_fract: " + cld_fract);
	System.out.println("frveg: " + frveg);
	System.out.println("xlai: " + xlai);
	System.out.println("epsf: " + epsf);
	System.out.println("albedo_fflag: " + albedo_fflag);
	System.out.println("albf: " + albf);
	System.out.println("stmtype: " + stmtype);
	System.out.println("volrel: " + volrel);
	System.out.println("rmin: " + rmin);
	System.out.println("rcut: " + rcut);
	System.out.println("wilt: " + wilt);
	System.out.println("vegheight: " + vegheight);
	System.out.println("width: " + width);
	System.out.println("steady: " + steady);
	System.out.println("ci: " + ci);
	System.out.println("co: " + co);
	System.out.println("coz_sfc: " + coz_sfc);
	System.out.println("coz_air: " + coz_air);
	System.out.println("rks: " + rks);
	System.out.println("cosbyb: " + cosbyb);
	System.out.println("thmax: " + thmax);
	System.out.println("psis: " + psis);
	System.out.println("mintemp: " + mintemp);
	System.out.println("maxtemp: " + maxtemp);
	System.out.println("beta: " + beta);
	System.out.println("b1: " + b1);
	System.out.println("b2: " + b2);
	System.out.println("psice: " + psice);
	System.out.println("sc: " + sc);
	System.out.println("zp: " + zp);
	System.out.println("frhgt: " + frhgt);
	System.out.println("frzp: " + frzp);
	System.out.println("rkocap: " + rkocap);
	System.out.println("rccap: " + rccap);
	System.out.println("rzcap: " + rzcap);
	System.out.println("volini: " + volini);
	System.out.println("zstini: " + zstini);
	System.out.println("nobs_ptq: " + nobs_ptq);
	System.out.println("nobs_wind: " + nobs_wind);
	System.out.println("station_height: " + station_height);
	System.out.println("ugs: " + ugs);
	System.out.println("vgs: " + vgs);*/
// java.io.FileOutputStream fout = new java.io.FileOutputStream("par.out");
// java.io.PrintStream p = new java.io.PrintStream(fout);
// java.io.FileOutputStream fout2 = new java.io.FileOutputStream("wind.out");
// java.io.PrintStream p2 = new java.io.PrintStream(fout2);
        // DEBUG print sounding data input arrays
/*
		
  for (int q=0; q<ps.length; q++) { p.println("ps[" + q + "] " + ps[q]); }
  for (int q=0; q<ts.length; q++) { System.out.println("ts[" + q + "] " + ts[q]); }
  for (int q=0; q<dep.length; q++) { System.out.println("dep[" + q + "] " + dep[q]); }
  for (int q=0; q<dd0.length; q++) { System.out.println("dd0[" + q + "] " + dd0[q]); }
  for (int q=0; q<ff0.length; q++) { System.out.println("ff0[" + q + "] " + ff0[q]); }
  for (int q=0; q<zh.length; q++) { System.out.println("zh[" + q + "] " + zh[q]); }
		  
		/* Type declaration of variables */

 /*
		 * Note that they are *all* initialised, which probably isn't necessary,
		 * but it's being done as a precaution; fortran assigns a new (numeric)
		 * variable the value of zero, in C it is undefined.
         */
        double advgt = 0, aepsi = 0, ahum = 0, albdoe = 0, angl = 0,
                aptemp = 0, atemp = 0, awind = 0, bulk = 0, caf = 0,
                capaci = 0, ccan = 0, cf = 0, cg = 0, cha = 0, chf = 0,
                chg = 0, chgt = 0, clktam = 0, clktpm = 0, cp = 0,
                ctheta = 0, del = 0, delt = 0, delta = 0, deltaz = 0,
                deltvst = 0, dhet = 0, dqdt = 0, dtx = 0, dty = 0,
                dzeta = 0, embar = 0, evap = 0, fc = 0, fco2 = 0,
                fglobal = 0, flux_plant = 0, fluxgd = 0, fpsie = 0,
                frco2 = 0, fs = 0, fst = 0, gam = 0, grav = 0, heat = 0,
                het = 0, hf = 0, hg = 0, hgt = 0;
        double[] abstbl = new double[51];
        double[] bsctbl = new double[51];
        double[] delz = new double[9];
        double[] dqdt2 = new double[51];
        double[] gm = new double[51];
        int ifirst = 0, jcap = 0, kqflag = 0, nlvls = 0, ntrp = 0;
        double o_pot_tmp = 0, oshum = 0, otemp = 0, ps1 = 0, psicm = 0,
                psie = 0, psig = 0, psim = 0, psist = 0, psisup = 0,
                psiwc = 0, psix = 0, ptime = 0, ptm100 = 0, qa = 0, qaf = 0,
                qstf = 0, qstg = 0, r = 0, rad = 0, raf = 0, realtm = 0,
                resist = 0, rkw = 0, rl = 0, rlelf = 0, rlf = 0, rlg = 0,
                rlpsi = 0, rnet = 0, rnetf = 0, rnetg = 0, rs = 0,
                rscrit = 0, rsf = 0, rsg = 0, rst = 0, rtranw = 0,
                rzascr = 0, satam = 0, satpm = 0, sigf = 0, sigma = 0,
                sol = 0, sum = 0, sumo3 = 0, sumw = 0, swave = 0, ta = 0,
                taf = 0, tdif_50 = 0, tdif_s = 0, tf = 0, tg = 0, thv = 0,
                tscren = 0, tstar = 0, uaf = 0, ustar = 0, uten = 0,
                vfl = 0, volist = 0, volrmv = 0, w2g = 0, wgg = 0, wpsi = 0,
                xlef = 0, xleg = 0, xmax = 0, xmod = 0, xpdif = 0,
                xtdif = 0, y = 0, ypdif = 0, ytdif = 0, za = 0, zcount = 0,
                zg = 0, zst = 0;
        double[] q_fine = new double[52];
        double[] qd = new double[51];
        double[] qn = new double[52];
        double[] qq = new double[51];
        double[] scatbl = new double[51];
        double[] t = new double[51];
        double[] t_fine = new double[52];
        double[] td = new double[51];
        double[] tt = new double[10];
        double[] u = new double[51];
        double[] u_fine = new double[52];
        double[] ud = new double[51];
        double[] ug = new double[51];
        double[] ugd = new double[51];
        double[] v = new double[51];
        double[] v_fine = new double[52];
        double[] vd = new double[51];
        double[] vg = new double[51];
        double[] vgd = new double[51];
        double[] xfun = new double[10];
        double[] z = new double[10];
        double[] zi = new double[51];
        double[] zk = new double[51];

        /* Non-standard (for fortran) naming */
        double[] km = new double[51];
        double lambda = 0, kappa = 0, le = 0, karman = 0, lwdn = 0, mol = 0;

        /* Some sort-of localish variables */
 /* advect */
        double dtdx = 0, dtdy = 0;
        int advect_pass = 0;

        /* air */
        double cdelt = 0, deltx = 0, tdel = 0, ttop = 0, zmix = 0;

        /* albedo */
        int algflg = 0, albedo_pass = 0, albflg = 0, alfflg = 0;

        /* avr */
        double[] av_array = new double[5];
        int average_pass = 0;
        double average_sum = 0;

        /* ball */
        char limit;
        double akco = 0, akoo = 0, oxy = 0, par = 0, aaa = 0, aj = 0,
                ajmax = 0, ajpar = 0, akc = 0, ako = 0, alpi = 0, ampar = 0,
                bbb = 0, ccc = 0, ccomp = 0, cii = 0, css = 0, gs = 0,
                parq = 0, qbnd = 0, resp = 0, respo = 0, rhl = 0, tau = 0,
                tauo = 0, vc = 0, vcmax = 0, vcmaxo = 0, we = 0, wr = 0,
                ws = 0;
        int ball_repeat = 0;

        /* below */
        double dummy = 0, term1 = 0, term2 = 0, term3 = 0;
        double[] dtdt = new double[9];
        double[] te = new double[10];
        double[] ttt = new double[10];

        /* bri */
        double cr1 = 0, mol1 = 0, phih = 0, psimnew = 0, psimzten = 0,
                radcor = 0, t2 = 0, tdif = 0, tsurf = 0, wg1 = 0,
                zoverl = 0, ztenoverl = 0;

        /* calc */
        double lat_in_radians = 0, out_time_intv = 0;

        /* capac */
        double addit = 0, capini = 0, caprat = 0, fpm = 0, voliso = 0,
                volrmo = 0, zaaa = 0, zuh = 0, zadd = 0;
        int idel = 0, itrap = 0, ixcap = 0, jdel = 0, ncap = 0;

        /* co2flx */
        double pex = 0, rafcanopy = 0, rroe = 0;

        /* flux */
        double evap_smooth = 0, mixedevapflux = 0, mixedradiotemp = 0,
                vegnevapflux = 0;

        /* Gettbl */
        double[] gt_a = new double[47];
        double[] gt_b = new double[47];
        double[] gt_c = new double[47];

        /* gtemp */
        double a1 = 0, deltax = 0, engbal = 0, xdif = 0;

        /* hot */
        double mixedheatflux = 0, vegnheatflux = 0;

        /* Input */
        double bscat = 0, bscatd = 0, dad = 0, de = 0, df = 0, gmt = 0,
                hi = 0, radian = 0, rlpath = 0, rm = 0, s = 0, sdec = 0,
                sheat = 0, siga = 0, skonst = 0, solel = 0, tabs = 0,
                tabsd = 0, tscat = 0, tscatd = 0, xser = 0;
        double[] ss = new double[13];
        int imol = 0, input_pass = 0, jmo = 0, imo1 = 0;
        int[] md = new int[13];

        /* intpol */
        double fix = 0, tix = 0, utdz = 0, vtdz = 0, zz = 0, zzday = 0;
        int intpol_pass = 0;

        /* main */
        double tmod = 0;

        /* mom */
        double abdu = 0, abdv = 0, bx1 = 0, krad = 0, okmax = 0, rc = 0,
                sb = 0, uif2 = 0;
        double[] bx = new double[47];
        double[] cr = new double[47];
        double[] dtw = new double[47];
        double[] duw = new double[47];
        double[] dvw = new double[47];
        double[] ri = new double[47];
        double[] wgeos = new double[47];
        double[] ok = new double[47];
        int imax = 0, imax1 = 0, mom_pass = 0;

        /* momday */
        double deltzi = 0, du50 = 0, dudt50 = 0, dv50 = 0, dvdt50 = 0,
                over = 0, smf = 0, smfh = 0, xxx = 0;
        double[] dudt = new double[51];
        double[] dvdt = new double[51];
        int lnugtjooclwwfli = 0, momday_pass = 0;

        /* netrad */
        double lwup = 0, vegnshortwave = 0;

        /* output */
        double air_leaf_t = 0, ccan_concn = 0, co2_flux = 0, g_flux = 0,
                stom_r = 0, undefined = 0, water_use_eff = 0;
        int i_columns = 0, output_pass = 0;

        /* ozone */
        double fbare = 0, fg = 0, fleaf = 0, fmeso = 0, rag = 0, rroz = 0,
                rtot = 0, the_time = 0;

        /* pslcal */
        double hbar = 0, q_boundary = 0, units1 = 0, units2 = 0, zeff = 0,
                ztot = 0;
        int pslcal_pass = 0, init_pslcal = 0;

        /* psoil */
        double corr = 0, oma = 0, sffk = 0, xincr = 0, zeta = 0;
        int ii = 0, psoil_pass = 0;

        /* slope */
        double dipan = 0, dipaz = 0;

        /* snding */
        double dydx = 0, dydxn = 0, eqt = 0, height_at_ntrp = 0,
                pot_50 = 0, precip_h20 = 0, pres_50 = 0, sum_precip_h20 = 0,
                tbar = 0, tdew = 0, vert_spacing = 0;
        double[] derivs;
        double[] ew = new double[51];
        double[] gmq = new double[51];
        double[] pot_temp = new double[51];
        double[] qs = new double[51];
        double[] ucomp = new double[51];
        double[] vcomp = new double[51];

        /* start */
        char cclass;
        char[] soiltype = new char[MAX_STRING_LENGTH];
        char[] planttype = new char[MAX_STRING_LENGTH];
        int index_soils = 0, index_veggies = 0, num_of_veggies = 0,
                num_soils = 0, ob_flag = 0, zo_flag = 0;

        /* vegflx */
        double rprime = 0;

        /* veghot */
        double aveg = 0, hfn = 0;

        /* vegrad */
        double tzero4 = 0, tzero5 = 0;

        /* vegvel */
        double cdl = 0, rhg = 0, rhocp = 0, sdl = 0, sigalf = 0, qsta = 0,
                rmratiodif = 0, tfn = 0, xlefn = 0;
        int vegvel_pass = 0;

        /* vel */
        double ch_obst_hgt = 0, chi20 = 0, chia = 0, chio = 0, cmh = 0,
                cmw = 0, fm = 0, fobst_hgt = 0, ft_obst_hgt = 0, ft20 = 0,
                ften = 0, kmm = 0, ks = 0, kw = 0, kx = 0, ptmp20 = 0,
                reflev = 0, rekust = 0, robst_patch = 0, rtrano3 = 0,
                rtrans = 0, rtrans_patch = 0, rtransw_patch = 0,
                rza_obst = 0, rza_obst_hgt = 0, rzazo = 0, sa = 0,
                so10m = 0, sobst_hgt = 0, srfczo = 0, sten = 0, t_ft = 0,
                u_patch = 0, uscrn = 0, ustar_patch = 0;
        int init_vel = 0, vel_pass = 0;

        /* water */
        double c11 = 0, c22 = 0, c33 = 0, c44 = 0, evai = 0, evas = 0,
                evax = 0, per = 0, win = 0, ww1 = 0, ww2 = 0, ww3 = 0;

        /* Other crap used in more than one routine */
        double actsqrt = 0, aroot = 0, b1_p = 0, b2_p = 0,
                bareevapflux = 0, bareheatflux = 0, barenetradn = 0,
                bareradiotemp = 0, bowen = 0, broot = 0, chax = 0,
                croot = 0, effdec = 0, fbscat = 0, fract = 0, fract2 = 0,
                ftabs = 0, ftscat = 0, h = 0, height = 0, hrangl = 0,
                mixednetradn = 0, oldtmp = 0, path = 0, pes = 0, pot_s = 0,
                psihnew = 0, rair = 0, rrtot = 0, rstdiv = 0, sig = 0,
                sinsolelslope = 0, slb = 0, solelslope = 0, solsin = 0,
                t1 = 0, unscaled_raf = 0, ustar1 = 0, vegnnetradn = 0,
                vegnradiotemp = 0, vfl2 = 0, x = 0, x1 = 0, zten = 0,
                thick_mx = 0, thick = 0, dz = 0, dz_advect = 0, uif = 0,
                uif_mom = 0, du = 0, dv = 0, sgma = 0, sgma_capac = 0,
                rha = 0, a_sh = 0, c = 0, b = 0;
        double[] wind = new double[47];
        double[] du_momday = new double[51];
        double[] dv_momday = new double[51];
        double[] a_gtemp = new double[6];
        double[] b_mom = new double[47];
        int dual_regime = 0, incr = 0, j = 0, k = 0, l = 0, n = 0,
                nlvl1 = 0, i_header_output = 0;

        /* New indices */
        double wbpt = 0;
        double dp_average = 0;
        double t_average = 0;
        double tat850, tat700, tat500, tdat850, tdat700, wsat850, wsat500, swcorr;
        double uat850, vat850, uat500, vat500, wdir500, wdir850;
        double ps_show, p850, thtlcl, guess, pr_tmst, p500;

        /*update by Yannis Konstas(start)*/
        double ef1 = 0, ef2 = 0, ed1 = 0, edN = 0, edJ = 0, edA = 0, edB = 0, itwu = 0, Ws = 0, Nn = 0;

        //calculations for Ed
        edA = 12.0d - 0.0569d * xlat - 0.000202d * (xlat * xlat)
                + 0.00000825d * (xlat * xlat * xlat) - 0.000000315d * (xlat * xlat * xlat * xlat);

        edB = 0.123d * xlat - 0.000310d * (xlat * xlat) + 0.0000008d * (xlat * xlat * xlat)
                + 0.000000499d * (xlat * xlat * xlat * xlat);

        edN = 0.945 * (edA + edB * Math.pow(Math.sin(Math.PI * ((double) ((double) dayOfTheYear + 10.0d) / 365.0d)), 2.0d));
        edN *= 1.0d / 24.0d;

//                System.out.println("day = "+ dayOfTheYear+ " xlat = "+ xlat+" edA = "+edA +" edB = "+edB + " edN = "+edN);
        /*update by Yannis Konstas(end)*/
 /* File pointers */
        // FILE *o_model, *LUT;

        /* More variables declared in main */
        char[] fractn = new char[MAX_STRING_LENGTH];
        int stabcriteria = 0, init = 0, no_rows = 0, ionce = 0,
                monce = 0;
        double zo_patch = 0, old_ahum = 0, ycount = 0, time = 0;
        double[] zls = new double[51];
        int out_count = 0;

        /* These are mine: dajr */
        int file_error;
        int i;

        /* Just to be extra safe, initialise all the arrays to zero. */
        for (i = 0; i < 5; i++) {
            av_array[i] = 0;
        }

        for (i = 0; i < 6; i++) {
            a_gtemp[i] = 0;
        }

        for (i = 0; i < 9; i++) {
            delz[i] = 0;
            dtdt[i] = 0;
        }

        for (i = 0; i < 10; i++) {
            tt[i] = 0;
            xfun[i] = 0;
            z[i] = 0;
            te[i] = 0;
            ttt[i] = 0;
        }

        for (i = 0; i < 13; i++) {
            md[i] = 0;
            ss[i] = 0;
        }

        for (i = 0; i < 47; i++) {
            gt_a[i] = 0;
            gt_b[i] = 0;
            gt_c[i] = 0;
            bx[i] = 0;
            cr[i] = 0;
            dtw[i] = 0;
            duw[i] = 0;
            dvw[i] = 0;
            ri[i] = 0;
            wgeos[i] = 0;
            ok[i] = 0;
            wind[i] = 0;
            b_mom[i] = 0;
        }

        for (i = 0; i < 51; i++) {
            abstbl[i] = 0;
            bsctbl[i] = 0;
            dqdt2[i] = 0;
            gm[i] = 0;
            qd[i] = 0;
            qq[i] = 0;
            scatbl[i] = 0;
            t[i] = 0;
            td[i] = 0;
            u[i] = 0;
            ud[i] = 0;
            ug[i] = 0;
            ugd[i] = 0;
            v[i] = 0;
            vd[i] = 0;
            vg[i] = 0;
            vgd[i] = 0;
            zi[i] = 0;
            zk[i] = 0;
            zls[i] = 0;
            km[i] = 0;
            dudt[i] = 0;
            dvdt[i] = 0;
            ew[i] = 0;
            gmq[i] = 0;
            pot_temp[i] = 0;
            qs[i] = 0;
            ucomp[i] = 0;
            vcomp[i] = 0;
            du_momday[i] = 0;
            dv_momday[i] = 0;
        }

        for (i = 0; i < 52; i++) {
            q_fine[i] = 0;
            qn[i] = 0;
            t_fine[i] = 0;
            u_fine[i] = 0;
            v_fine[i] = 0;
        }

        // System.out.println("Variables defined");

        /* Define the "global" variables here */
        ionce = 0;
        monce = 0;
        ycount = 0;
        time = 0;
        stabcriteria = 0;
        init = 1;

        /* Set up some variables, done in "block" */
        y = 1.0;
        xmod = 0.0;
        sigf = 0.0;
        hg = 0.0;
        ahum = 0.0;
        rnet = 0.0;

        for (i = 1; i <= 21; i++) {
            qd[i] = 0;
        }

        chgt = 0.0;
        ustar = 0.0;
        tstar = 0.0;
        heat = 0.0;
        hgt = 50.0;
        za = 50.0;
        delt = 1;
        ctheta = 1;
        dhet = 0.0;
        evap = 0.0;
        sigma = 5.6521e-8;
        le = 2.5e6;
        karman = 0.4;
        grav = 9.78;
        r = 287.5;
        rad = 1.6e-5;
        cp = 1004.832;
        delta = 90;
        mol = 0;
        bulk = 0;
        ifirst = 0;
        nlvls = 5;
        jcap = 1;

        /* Open the output file */
        // o_model=fopen("o_model.dat","w");
        // if(o_model==null)
        // {
        // System.out.println("Error opening output file\n");
        // throw new Exception();
        // }

        /* Read and check data */
        /**
         * ******************************************
         * start ******************************************
         */

        /* Read in the control variables from the input file. */
 /* This is where the file reads were originally */
        zo_flag = 0;
        ob_flag = 0;
        cclass = 'u';

        /* Include the vegetation and soils databases in the in calculations */
        if (xlai == 0) {
            xlai = 1;
        }

        xlef = 1;

        /* This is where the soil and custom plant parameters were read. */
 /* Make some more conversions from one unit to another */
        btemp = (btemp + 273.15);

        // TCS Why outtt is multiplied by 60, and then divided by 60 the next time it is used,
        // is beyond my capacity for judging.
        //outtt = outtt * 60;
        co = (co * 1.0e-6);
        ci = (ci * 1e-6);
        mintemp = mintemp + 273;
        maxtemp = maxtemp + 273;

        /* The relevant data for rough.for */
        if (frveg == 0.0) {
            cclass = 'B';
        }

        if (frveg == 100.0) {
            cclass = 'V';
        }

        if ((frveg < 100.0) || (frveg > 0.0)) {
            cclass = 'P';
        }

        if (zo != 0) {
            zo_flag = 1;
        }

        if (obst_hgt != 0) {
            ob_flag = 1;
        }

        /* Finished reading the data */
        // System.out.println("Variables initialized - start routines finished");
        /**
         * ******************************************
         * end of start ******************************************
         */
        /**
         * ******************************************
         * snding ********************************************? Read sounding
         * and call spline to interpolate
         */

        /* Data vertical spacing */
        vert_spacing = 250;
        /* 250 meter intervals */
        deltaz = vert_spacing;

        /* This is where the sounding data was read in */
 /* Do some initial calculations on the soundings. */
        height = 0;

        for (j = 1; j <= nobs_ptq; j++) {
            i = j - 1;
            tdew = ts[j] - dep[j] + C_TO_K;
            /* Dew point in kelvin */
            ew[j] = (6.11 * Math.exp((2.49e6 / 461.51)
                    * (1 / C_TO_K - 1 / tdew)));
            /* es in millibars */
            qs[j] = (0.622 * ew[j] / ps[j]);
            /* Specific humidity in g/kg */

            pot_temp[j] = ((ts[j] + C_TO_K) * Math.pow((1000.0 / ps[j]), 0.286));
            /* Theta in kelvin */

            if (j > 1) {
                tbar = ((ts[i] + ts[i + 1]) / 2) + C_TO_K;
                /* Average temperature in kelvin */
                thick = (287 * (tbar / 9.8) * Math.log(ps[i] / ps[i + 1]));
                /* Thickness in meters */
                height = height + thick;
                /* Height of pressure level above station */
                zls[i + 1] = height;
                gm[i] = (pot_temp[i + 1] - pot_temp[i]) / thick;
                /* d(theta)/dZ */
                gmq[i] = (qs[i] - qs[i + 1]) / thick;
                /* d(q)/dZ */
                precip_h20 = (-0.622 / GRAV * 10.0
                        * ((ew[j - 1] * ps[j] - ew[j] * ps[j - 1])
                        * Math.log(ps[j] / ps[j - 1])
                        / (ps[j] - ps[j - 1]) + ew[j] - ew[j - 1]));
                sum_precip_h20 = sum_precip_h20 + precip_h20;
                omega = sum_precip_h20;
            }
        }

        /* Winds - note the way that they're recorded by the weather service */
        zh[1] = 0;

        for (i = 2; i <= nobs_wind; i++) {
            zh[i] = ((zh[i] - station_height) * 1000) / FT_TO_METRES;
        }

        /* Calculate u and v velocity components */
        for (j = 1; j <= nobs_wind; j++) {
            ucomp[j] = ((-ff0[j] * Math.cos((90.0 - dd0[j]) / RADIANS_TO_DEGREES))
                    / KTS_TO_METRES);
            vcomp[j] = ((-ff0[j] * Math.sin((90.0 - dd0[j]) / RADIANS_TO_DEGREES))
                    / KTS_TO_METRES);
        }

        /* Numpber of levels it's possible to interpolate (ntrp) are */
        height_at_ntrp = _min(height, zh[nobs_wind]);
        ntrp = (int) ((height_at_ntrp - 50.0) / vert_spacing);

        /* Interpolate at 'h' height intervals using cubic splines. */
 /* Numerical recipes in C */
 /* 1. Temperature */

 /*
		 * Allocate space for derivatives. Note extra one, thanks to unit-offset
		 * indexing. Fortunately (or rather not) "numerical recipes in C" (from which
		 * the spline interpolation is ripped) insists on using unit-offsets, despite
		 * the fact that no-one in the whole world of C has ever done so before or
		 * since.
         */
        // derivs=calloc((nobs_ptq+1),sizeof);
        derivs = new double[nobs_ptq + 1];

        if (derivs == null) {

            // System.out.println("Error allocating memory for derivatives\n");
            throw new Exception("Error allocating memory for derivatives\n");
        }

        /* Calculate derivatives at the boundaries */
        dydx = (pot_temp[2] - pot_temp[1]) / (zls[2] - zls[1]);
        dydxn = (pot_temp[nobs_ptq] - pot_temp[nobs_ptq - 1])
                / (zls[nobs_ptq] - zls[nobs_ptq - 1]);

        /* Calculate array of derivatives */
        spline(zls, pot_temp, nobs_ptq, dydx, dydxn, derivs);

        /* Call splint to get actual value at 50m and subsequent N metre intervals */
        for (i = 0; i <= ntrp; i++) {
            h = 50 + (vert_spacing * i);
            zi[i + 1] = h;

            // splint(zls, pot_temp, derivs, nobs_ptq, h, &td[i+1]);
            td[i + 1] = splint(zls, pot_temp, derivs, nobs_ptq, h);
        }

        /* Free up the space taken by the list of derivatives */
        // free(derivs);
        derivs = null;

        // System.gc();
        // System.out.println("Temperature soundings calculated");

        /* 2. Moisture */
 /* Why on earth his is allocated again, I have no idea. */

 /*
		 * Why on earth this next section is not incorportated into the previous
		 * one, I have no idea either.
         */

 /*
		 * Allocated memory for derivatives, again, even though there's the same
		 * number as before
         */
        // derivs=calloc((nobs_ptq+1),sizeof);
        derivs = new double[nobs_ptq + 1];

        if (derivs == null) {

            // System.out.println("Error allocating memory for derivatives\n");
            throw new Exception("Error allocating memory for derivatives\n");
        }

        /* Calculate the derivatives at the boundaries */
        dydx = (qs[2] - qs[1]) / (zls[2] - zls[1]);
        dydxn = (qs[nobs_ptq] - qs[nobs_ptq - 1])
                / (zls[nobs_ptq] - zls[nobs_ptq - 1]);

        /* Calculate the array of derivatives. */
        spline(zls, qs, nobs_ptq, dydx, dydxn, derivs);

        /* Call splint to get actual value at 50m and subsequent N metre intervals */
        for (i = 0; i <= ntrp; i++) {
            h = 50 + (vert_spacing * i);
            qd[i + 1] = splint(zls, qs, derivs, nobs_ptq, h);
        }

        /* Free up the space taken by list of derivatives */
        // free(derivs);
        derivs = null;

        // System.gc();
        // System.out.println("Moisture soundings calculated");

        /* 3. Winds. */
 /* Calculate the u component */

 /*
		 * Create and array of floats for the derivatives. Well gosh we actually
		 * needed to reallocate the space this time. # of wind levels is not
		 * necessarily equal to the number of ptq levels
         */
        // derivs=calloc((nobs_wind+1),sizeof);
        derivs = new double[nobs_wind + 1];

        if (derivs == null) {

            // System.out.println("Error Allocating memory for derivatives\n");
            throw new Exception("Error Allocating memory for derivatives\n");
        }

        /* Calculate derivatives at the boundaries */
        dydx = (ucomp[2] - ucomp[1]) / (zh[2] - zh[1]);
        dydxn = (ucomp[nobs_wind] - ucomp[nobs_wind - 1])
                / (zh[nobs_wind] - zh[nobs_wind - 1]);

        /* Calculate array of derivatives */
        spline(zh, ucomp, nobs_wind, dydx, dydxn, derivs);

        /* Call splint to get the actual value at 50m and subsequent 250m intervals */
        for (i = 0; i <= ntrp; i++) {
            h = 50 + (vert_spacing * i);
            ud[i + 1] = splint(zh, ucomp, derivs, nobs_wind, h);
        }

        /*
		 * Gee, freeing up that array again, just to allocate at again
		 * with the same size
         */
        // free(derivs);
        derivs = null;

        // System.gc();

        /*
		 * Same as above, but for the v component. As before with the PTQ, this
		 * should ne incorporated into a single loop along with the u component.
         */
        // derivs=calloc((nobs_wind+1),sizeof);
        derivs = new double[nobs_wind + 1];

        if (derivs == null) {

            // System.out.println("Error allocating memory for derivatives\n");
            throw new Exception("Error allocating memory for derivatives\n");
        }

        /* Calculate derivates at the boundaries */
        dydx = (vcomp[2] - vcomp[1]) / (zh[2] - zh[1]);
        dydxn = (vcomp[nobs_wind] - vcomp[nobs_wind - 1])
                / (zh[nobs_wind] - zh[nobs_wind - 1]);

        /* Calculate array of derivatives */
        spline(zh, vcomp, nobs_wind, dydx, dydxn, derivs);

        /* Call splint to get the actual value at 50m and subsequent 250m intervals */
        for (i = 0; i <= ntrp; i++) {
            h = 50 + (vert_spacing * i);
            vd[i + 1] = splint(zh, vcomp, derivs, nobs_wind, h);
        }

        /* Free up that array again. But hey, we actually should this time. */
        // free(derivs);
        derivs = null;

        // System.gc();

        /*
		 * Calculate the lapse rate in the first layer and ew at screen level
		 * for use in the calc of screen level sat'n spec humidity
         */
        atemp = 50 * (ts[2] - ts[1]) / zls[2] + ts[1] + C_TO_K;
        tscren = ts[1] + C_TO_K;
        /* Screen temperature */

 /* Assign value to all elements of ew[] */
        for (i = 0; i <= 50; i++) {
            ew[i] = (6.11 * Math.exp((2.49e6 / 461.51)
                    * (1 / C_TO_K - 1 / tscren)));
        }

        oshum = (0.622 * ew[1] / ps[1]);
        ahum = qs[1];
        old_ahum = qs[1];
        ps1 = ps[1];

        /* Changes 2/10/92 */
        pres_50 = (ps[1] * Math.exp(-9.8 * 50 / (287 * (ts[1] + C_TO_K))));
        pot_50 = (atemp * Math.pow((1000 / pres_50), 0.286));
        aptemp = pot_50;
        o_pot_tmp = pot_50;
        tdif_50 = pot_50 - atemp;
        tdif_s = (tdif_50 - 0.5);

        /*
		 * Here's some left-over crap from goodness only knows when. I'm leaving it
		 * in as comments in case someone might want to refer to it at some point.
		 * (This is left as fortran)
         */

 /*
		 * tdeww=tdew=273.15
		 * expt=7.5*tdeww/(273.3+tdeww)
		 * eww=6.11*10**expt
		 * ewww=qs(j)*ps(j)/(0.378*qs(j)+0.622)
		 * 
		 * rhoa=ps(j)*100./(287.*tbar) ! Kell's precip water calc
		 * asum=asum+(qs(i)+qs(i+1))/2.*.001*rhoa*thick*100
         */
        /**
         * ******************************************
         * The end of snding, finally.
         * ******************************************
         */
        // System.out.println("Winds calculated. Done with soundings");

        /* Some basic calculations */
        /**
         * ******************************************
         * calc ******************************************
         */
        k = (int) xlat;
        xlat = ((xlat - k) / 0.6 + k);
        k = (int) xlong;
        xlong = ((xlong - k) / 0.6 + k);

        /* Convert lat to radians, calculate coliolis force */
        lat_in_radians = xlat * DEGREES_TO_RADIANS;
        cf = 2 * ROT_RATE_EARTH * Math.sin(lat_in_radians);

        /* Time conversions */
 /* TCS All times are passed in as total seconds, so no conversions are needed */
        // timend=_dectim(timend); /* Convert to decimal time */
        // strtim=_dectim(strtim);	
        //out_time_intv = outtt / 60;	
        //out_time_intv = _dectim(out_time_intv); 
        out_time_intv = outtt;
        no_rows = (((int) timend - (int) strtim) / (int) out_time_intv) + 1;

        /**
         * Unused TCS satam=_dectim(clktam); satpm=_dectim(clktpm);
         */

        /* Store and initialise temperatures */
        oldtmp = atemp;
        t[1] = (atemp + 0.5);
        otemp = tscren - 2;
        frveg = frveg / 100;
        /* Vegetation %age as a fraction */

        /**
         * ******************************************
         * end of calc ******************************************
         */

        // System.out.println("Done basic calcs.");

        /* Set a geostrophic wind at the surface */
        /**
         * ******************************************
         * prfile ******************************************
         */

        /*
		 * Generate the daytime and nighttime vertical wind profiles of the
		 * geostrophic wind components at intervals of 250m from 50m (top of the
		 * surface layer) from the surface geostrophic winds.
		 * They can be calculated in three different ways; which is specified by
		 * the user.
         */
        j = 1;
        u_fine[1] = ud[1];
        v_fine[1] = vd[1];
        t_fine[1] = td[1];
        q_fine[1] = qd[1];

        for (i = 1; i <= 10; i++) {
            for (incr = 1; incr <= 5; incr++) {
                j = j + 1;
                u_fine[j] = _rlinear(ud[i + 1], ud[i], incr);
                v_fine[j] = _rlinear(vd[i + 1], vd[i], incr);
                t_fine[j] = _rlinear(td[i + 1], td[i], incr);
                q_fine[j] = _rlinear(qd[i + 1], qd[i], incr);
            }
        }

        /**
         * ******************************************
         * intpol ******************************************
         */

        /*
		 * Set a couple of variables the first time through. This is a really
		 * ugly way of doing it, but it's there in the fortran.
         */
        if (intpol_pass == 0) {
            zz = 50;
            zzday = 250;
            intpol_pass = 1;
        }

        vgd[1] = vgs;
        ugd[1] = ugs;

        /* Set the geostrophic wind component equal winds at 1050m and above */
        for (j = 5; j <= ntrp; j++) {
            vgd[j] = vd[j];
            ugd[j] = ud[j];
        }

        /*
		 * Make the winds above ntrp levels equal that at the ntrp level.
		 * Trap so that the winds are not zeroed above ntrp levels
         */
        for (i = ntrp; i <= 20; i++) {
            vgd[i] = vgd[ntrp];
            ugd[i] = ugd[ntrp];
        }

        /*
		 * Set up a weighting factor to interpolate from 50 to 1050m.
		 * Calculate vertical derivatives for u & v
         */
        fix = 4 * zzday + 50;
        vtdz = (vd[5] - vgs) / fix;
        utdz = (ud[5] - ugs) / fix;

        /* Calculate the geostrophic wind components in intermediate levels */
        for (i = 1; i <= 4; i++) {
            tix = (i - 1) * zzday + 50;
            ugd[i] = ugs + utdz * tix;
            vgd[i] = vgs + vtdz * tix;
        }

        /* Vertical derivative for the night-time */
        utdz = (ud[3] - ugs) / 550;
        vtdz = (vd[3] - vgs) / 550;

        /* Follow procedure again for the night-time (50-500m) */
        for (i = 1; i <= 20; i++) {
            tix = i * zz;
            ug[i] = ugs + utdz * tix;
            vg[i] = vgs + vtdz * tix;
        }

        for (i = 21; i <= 46; i++) {
            ug[i] = u_fine[i];
            vg[i] = v_fine[i];
        }

        /**
         * ******************************************
         * end of intpol ******************************************
         */
        // System.out.println("Done with interpolations");
        /**
         * ******************************************
         * advect ******************************************
         */

        /*
		 * Calculates the geostrophic temperature advection based on the thermal
		 * wind equation and the vertical distribution of geostrophic wind.
         */
 /* Another ugly initialization */
        if (advect_pass == 0) {
            dz_advect = 1000;
            advect_pass = 1;
        }

        /*
		 * The depth of the layer over which we calculate the geostrophic
		 * temperature advection is 1000m, z[5]-z[1].
         */
        dtdx = cf * otemp / (GRAV * dz_advect) * (vgd[5] - vgd[1]);
        dtdy = (-cf) * otemp / (GRAV * dz_advect) * (ugd[5] - ugd[1]);
        advgt = (-1) * (ugd[3] * dtdx + vgd[3] * dtdy);

        /*
		 * Assume that the actual temperature change is one half that of
		 * geostrophic temperature advection.
         */
        advgt = advgt / 2;

        /**
         * ******************************************
         * end of advect ******************************************
         */

        /* Calculate the initial value of the wind at 50m... vel, bri */
        awind = Math.sqrt(ud[1] * ud[1] + vd[1] * vd[1]);

        /**
         * ******************************************
         * end of prfile ******************************************
         */
        // System.out.println("Done with profile");

        /* The lookup table for transm */
        /**
         * ******************************************
         * gettbl ******************************************
         */

        /*
		 * This section of the code takes the precipitable water (omega) and
		 * calculates the transmission coefficient for absorption using a lookup
		 * table. Interpolate linearly on omega. The lookup table contains 46-entry
		 * tables for omega from 0 to 5 in increments of 0.5.
		 * 
		 * Note that the maximum value allowed for omega is 5.0.
		 * 
		 * Code also copies the scattering (scatbl) and the backscattering (bsctbl)
		 * tables from the file into the common block.
         */

 /*
		 * This part might be kind of message to convert into Java, since (I
		 * assume) there's no easy way to open a file. The alternative might be to
		 * hard-code the lookup-tables into arrays and get them that way.
         */
        // LUT=fopen("LUT.dat","r");

        /*
		 * if(LUT==null)
		 * {
		 * System.out.println("Error opening lookup table\n");
		 * throw new Exception("Error opening lookup table\n");
		 * }
         */
        i = (int) (1 + omega * 2);
        j = i - 1;

        if (i > 11) {
            i = 11;
        }

        int lineNum = 1;	// The "line number" of the lut file

        /*
		 * For no precipitable water use the 1st 46 values, otherwise skip j
		 * blocks of 46 table entries and read appropriate values
         */
 /* Why the reading is duplicated, I have no idea, */
        if (j == 0) {
            for (l = 1; l <= 46; l++) {
                abstbl[l] = LUT[lineNum][0];
                scatbl[l] = LUT[lineNum][1];
                bsctbl[l] = LUT[lineNum][2];
                lineNum++;

                // IOIOIO //
                // file_error=fscanf(LUT,"%f %f %f",&abstbl[l],&scatbl[l],&bsctbl[l]);
                // if(file_error!=3)
                // {
                // System.out.println("Error reading lookup table\n");
                // throw new Exception();
                // }
            }
        } else {

            for (k = 1; k <= j; k++) {
                for (l = 1; l <= 46; l++) {
                    gt_a[l] = LUT[lineNum][0];
                    gt_b[l] = LUT[lineNum][1];
                    gt_c[l] = LUT[lineNum][2];
                    lineNum++;
                }
            }

            for (l = 1; l <= 46; l++) {
                abstbl[l] = LUT[lineNum][0];
                scatbl[l] = LUT[lineNum][1];
                bsctbl[l] = LUT[lineNum][2];
                lineNum++;

                // IOIOIO //
                // file_error=fscanf(LUT,"%f %f %f",&abstbl[l],&scatbl[l],&bsctbl[l]);
                // if(file_error!=3)
                // {
                // System.out.println("Error reading lookup table\n");
                // throw new Exception();
                // }
            }
        }

        /*
		 * For a maximum omega no interpolation is performed but for all other
		 * values of omega interpolation is done using weighted factors fract2 and
		 * fract. Values are written to channel 69.
         */

 /*
		 * Except nothing is written anywhere; there are just a bunch of loops
		 * with nothing (that isn't commented out) inside them. Still, here they all
		 * are, pointless, I know.
         */
        if (i == 11) {

            /*
                 * TCS - This just prints so it's commented out
                 * 
                 * printf("Lookup table for omega: %10.3f\n",omega);
                 * for(l=1;l<=46;l++)
                 * {
                 * printf("%12.7f %12.7f %12.7f\n",abstbl[l],scatbl[l],bsctbl[l]);
                 * }
             */
        } else {
            for (l = 1; l <= 46; l++) {
                gt_a[l] = LUT[lineNum][0];
                gt_b[k] = LUT[lineNum][1];
                gt_c[l] = LUT[lineNum][2];
                lineNum++;

                // IOIOIO //
                // file_error=fscanf(LUT,"%f %f %f",&gt_a[l],&gt_b[k],&gt_c[l]);
                // if(file_error!=3)
                // {
                // System.out.println("Error reading LUT.\n");
                // throw new Exception();
                // }
            }

            fract = 2 * omega - j;
            fract2 = 1 - fract;

            for (l = 1; l <= 46; l++) {
                abstbl[l] = abstbl[l] * fract2 + gt_a[l] * fract;
            }

            /*
                         * printf("GETTLB: Interpolation factors: %f %f\n",fract,fract2);
                         * printf("Lookup table for omega: %10.3f\n",omega);
             */

 /*
                         * TCS this does nothing
                         * for(l=1;l<=46;l++)
                         * {
                         * 
                         * printf("%17.7f %12.7f %12.7f\n",abstbl[l],scatbl[l],bsctbl[l]);
             */
            // }
            // END do nothing
        }

        // fclose(LUT);
        /**
         * ******************************************
         * end of gettbl ******************************************
         */
        // System.out.println("Done with lookup table");
        // System.out.println("---------- Starting diurnal loop");
        /**
         * * This is the start of the diurnal loop (timend-strtim) Nominally 22
         * hours *
         */

        /*
		 * Except it isn't really the start of the loop, the loop comes after
		 * the soil setup stuff.
         */
 /* Set up the soil */
        /**
         * ******************************************
         * psoil ******************************************
         */

        /*
		 * The usual ugly variable initialisation that probably isn't necessary,
		 * but it might be and I'd rather not take any chances.
         */
        if (psoil_pass == 0) {
            oma = 0.0000727;
            sffk = 4;
            psoil_pass = 1;
        }

        zeta = 0;
        dzeta = 1;

        /*
		 * Calc. lambda, kappa, and volumetric heat capacity of the ground (cg)
		 * See manual (HA!) for explanation of lambda and kappa.
         */

 /*
		 * Gillies, 1/11/95: The thermal inertia (tp) is entered in SI (Wm-2K-1)
		 * units. Convert to CGS to be able to use the regression equation derived
		 * from Sellers. This regression equation should be redone with MKS units.
         */
 /* Initiate dual ti option, tnc dajr marth 1996. */
        tp = (tp / 356.89466);
        /* Conversion to cal cm-1 K-1 s-1/2 */

        if (dual_ti == 'y') {
            tp = ti_a * f + ti_b;
        }

        /*
		 * More legacy stuff, commented out in the fortan, and still commented out
		 * here, but it's is here, because I never know what Toby might want to be able
		 * to see at some point in the future.
         */
 /* There are tp as a function of the surface moisture availability. Not used. */
 /* tp=0.04000*f+0.02000; /* Toby's generic soil. */
 /* tp=0.03000*f+0.03400; /* Toby's FIFE soil. */
 /* tp=0.04600*f+0.01400; /* Price's sand */
 /* tp=0.03890*f+0.01290; /* Price's clay. */
 /* tp=0.04241*f+0.01349; /* Price's AVG (snad, clay. */
 /* tp=0.03000*f+0.00690; /* Price's peat. */
 /* tp=0.04600*f+0.04500; /* Walnut gulch */
        lambda = (-0.00013 + 0.0502 * tp + 1.21 * tp * tp);
        kappa = (lambda * lambda) / (tp * tp);
        kappa = kappa / 10000;
        /* Now convert to MKS; tp entered in CGS. */
        lambda = (lambda * 418.68);
        cg = lambda / kappa;

        /*
		 * Find the best depth profile as a function of KAPPA.
		 * The following calcs are used in below. Creates a temp profile
		 * in the soil. Figure levels using scheme by Deardorff.
         */

 /*
		 * del is set to numerically stable value. Initial temps at substrate
		 * levels (nlvls) calculated based on a linear interpolation between
		 * otemp and btemp
         */
        nlvl1 = nlvls + 1;
        del = Math.sqrt(sffk * delta * kappa) / (Math.exp(dzeta) - 1);
        xincr = (btemp - otemp) / nlvls;
        tt[1] = otemp;
        tt[nlvl1] = btemp;

        for (ii = 2; ii <= nlvls; ii++) {
            tt[ii] = tt[ii - 1] + xincr;
        }

        /*
		 * Here the substrate spacing (depth) is being calculated based on the
		 * scale depth (h=1+z/d).
         */
        for (i = 1; i <= nlvl1; i++) {
            z[i] = (Math.exp(zeta) - 1) * del;
            xfun[i] = (1 + z[i] / del);
            zeta = zeta + dzeta;
        }

        /*
		 * Correction of lambda to account for reduction in temp wave in 1st
		 * soil layer. Consult manual for explanation.
         */
 /* There was an orphan "7678 continue" here in the fortran. */
        corr = 2 / (1 + Math.exp(-z[2] * Math.sqrt(oma / (2 * kappa))));
        lambda = corr * lambda;

        /* Initialise w2g and wgg, the substrate water budget parameters. */
        w2g = wmax * fsub;
        wgg = wmax * f;

        /**
         * ******************************************
         * end of psoil ******************************************
         */
        // System.out.println("Done with soil setup. Now I'm really starting the loop");

        /* Start of the main loop of the program */
        do {
            realtm = time + strtim;
            ptime = realtm / 3600.0;
            // System.out.println("** Loop iteration realtm " + realtm + " time " + time + " tmod " + tmod);

            // TCS I THINK THIS IS THE MODULUS FUNCTION, hope so : tmod=fmod(time, outtt);
            tmod = time % outtt;

            // System.out.println("** Loop iteration realtm " + realtm + " time " + time + " tmod " + tmod);
            /* Net radiation */
            /**
             * ******************************************
             * netrad *****************************************
             */

            /* Computes up and down longwave fluxes and the net radiation. */
            /**
             * ******************************************
             * input *****************************************
             */

            /*
                         * initialise some variables that may or may not be constant. Of
                         * course, they almost certainly are, but whatever.
             */
            if (input_pass == 0) {
                md[1] = 31;
                md[2] = 29;
                md[3] = 31;
                md[4] = 30;
                md[5] = 31;
                md[6] = 30;
                md[7] = 31;
                md[8] = 31;
                md[9] = 30;
                md[10] = 31;
                md[11] = 30;
                md[12] = 31;
                ss[1] = 0.967;
                ss[2] = 0.971;
                ss[3] = 0.982;
                ss[4] = 0.999;
                ss[5] = 1.016;
                ss[6] = 1.029;
                ss[7] = 1.034;
                ss[8] = 1.030;
                ss[9] = 1.018;
                ss[10] = 1.002;
                ss[11] = 0.985;
                ss[12] = 0.973;
                radian = 57.2957795;
                sdec = 0.39784988;
                siga = 279.9348;
                input_pass = 1;
            }

            skonst = (1.94 * 4.1868e4 / 60);

            /*
                         * y is counter for determining celestial mechanics for long, lat and
                         * time of year.
             */
            if (y == 1) {
                y = y + 1;
                imo1 = imo + 1;

                if (imo1 == 13) {
                    imo1 = 1;
                }

                /* Compute the solar distance factor */
                rm = (iday - 1) / (md[imo]);
                s = ss[imo] + rm * (ss[imo1] - ss[imo]);

                /* Set up the number of days in the year */
                k = iday;
                jmo = imo - 1;

                if (jmo >= 1) {
                    for (i = 1; i <= jmo; i++) {
                        k = k + md[i];
                    }
                }

                /* This wouldn't work for 1900 or 2100 */
 /* n=iyr/4; */
 /* n=iyr-n*4; */
 /* Amended to properly determine leap years */
                if (iyr % 100 == 0) {
                    n = iyr % 400;
                } else {
                    n = iyr % 4;
                }

                if ((n != 0) && (k >= 60)) {
                    k = k - 1;
                }

                dad = k - 1;

                /* Calculate the angular fraction of a year, convert to radians */
                df = (dad * 360 / 365.242);
                de = df / radian;

                /* Correction to declination caused by elliptical orbit. */
                sig = siga + df + 1.914827 * Math.sin(de)
                        - 0.079525 * Math.cos(de) + 0.019938 * Math.sin(de * 2)
                        - 0.00162 * Math.cos(de * 2);
                sig = sig / radian;

                /* Declination of the sun */
                effdec = Math.asin(sdec * Math.sin(sig));

                /* True solar noon */
                eqt = 0.12357 * Math.sin(de) - 0.004289 * Math.cos(de)
                        + 0.153809 * Math.sin(de * 2)
                        + 0.060783 * Math.cos(de * 2);
                eqt = eqt + 12;
            }

            /* Calculate time in greenwich mean time */
            gmt = ptime + tz;

            /* Calculate the solar hour angle in radians. */
            hrangl = 15 * (gmt - eqt) - xlong;
            hrangl = hrangl / radian;

            /* Now we can finally compute the solar elevation angle. */
            slb = xlat / radian;
            solsin = Math.sin(effdec) * Math.sin(slb)
                    + Math.cos(effdec) * Math.cos(slb) * Math.cos(hrangl);
            solel = Math.asin(solsin) * radian;

            if (slope > 0) {

                /**
                 * ******************************************
                 * sslope *****************************************
                 */

                /*
                                 * Calculates solar elevation angle and azimuth for sloping
                                 * terrain when slope is non-zero. Elevation for northwest corner
                                 * of square box is znw, etc. Grid size (called xmeshl) is in same
                                 * units as znw.
                 */

 /*
                                 * There was a whole bunch of stuff here just commented out, too
                                 * much for me to bother duplicating it.
                 */
                dipaz = (aspect / 57.2958);
                dipan = (slope / 57.2958);

                /* Computer solar elevation angle for sloping terrain. */
                sinsolelslope
                        = Math.cos(dipan)
                        * (Math.sin(slb) * Math.sin(effdec) + Math.cos(slb) * Math.cos(effdec) * Math.cos(hrangl))
                        + Math.sin(dipan)
                        * (Math.cos(dipaz)
                        * (Math.tan(slb) * (Math.sin(slb) * Math.sin(effdec) + Math.cos(slb) * Math.cos(effdec) * Math.cos(hrangl)) - Math.sin(effdec) / Math.cos(slb))
                        + Math.sin(dipaz) * Math.cos(effdec)
                        * Math.sin(hrangl));
                solelslope = Math.asin(sinsolelslope);

                if (solelslope <= 0.01) {
                    solelslope = 0.01;
                }

                /**
                 * ******************************************
                 * end of sslope *****************************************
                 */
                /**
                 * ******************************************
                 * albedo ******************************************
                 */

                /* Modified 4/26/01 to allow user to specify Albedoes */

 /*
                                 * If the albedoes are omitted in start then they're calculated,
                                 * otherwise the albedo of the ground and the foliage input are used to
                                 * calculate the weighted albedo (albdoe).
                                 * Note that albf depends on the solar angle (solsin); albg on wgg
                 */
                if (albedo_pass == 0) {
                    albflg = 0;
                    algflg = 0;
                    alfflg = 0;
                    albedo_pass = 1;
                }

                if (albflg == 0) {
                    if (albg == 0) {
                        algflg = 1;
                    }

                    if (albf == 0) {
                        alfflg = 1;
                    }

                    albflg = 1;
                }

                if (algflg == 1) {
                    albg = (0.25 - 0.20 * wgg / wmax);
                    /* Toby's value */

 /* albg=0.20-0.15*wgg/wmax; */
 /* Fudged dim */
 /* albg=0.30-0.20*wgg/wmax; */
 /* Fudged bright */
                }

                if (alfflg == 1) {

                    /* albf=0.025/(0.1+0.1*solsin)+(1-solsin*solsin)*0.1; */
                    albf = (0.032 / (0.1 + 0.1 * sinsolelslope)
                            + (1 - sinsolelslope * sinsolelslope) * 0.1);
                }

                sigf = 1 - Math.exp(-0.4 * xlai);
                albdoe = (sigf * albf + (1 - sigf) * albg) * frveg
                        + (1 - frveg) * albg;

                /**
                 * **********************************************
                 * end of albedo *********************************************
                 */
            } else {

                /* Albedo stuff all over again */
                if (albedo_pass == 0) {
                    albflg = 0;
                    algflg = 0;
                    alfflg = 0;
                    albedo_pass = 1;
                }

                if (albflg == 0) {
                    if (albg == 0) {
                        algflg = 1;
                    }

                    if (albf == 0) {
                        alfflg = 1;
                    }

                    albflg = 1;
                }

                if (algflg == 1) {
                    albg = (0.25 - 0.20 * wgg / wmax);
                    /* Toby's value */

 /* albg=0.20-0.15*wgg/wmax; */
 /* Fudged dim */
 /* albg=0.30-0.20*wgg/wmax; */
 /* Fudged bright */
                }

                if (alfflg == 1) {

                    /* albf=0.025/(0.1+0.1*solsin)+(1-solsin*solsin)*0.1; */
                    albf = (0.032 / (0.1 + 0.1 * solsin)
                            + (1 - solsin * solsin) * 0.1);
                }

                sigf = 1 - Math.exp(-0.4 * xlai);
                albdoe = (sigf * albf + (1 - sigf) * albg) * frveg
                        + (1 - frveg) * albg;

                /* End of repeated albedo stuff */
                sinsolelslope = solsin;
            }

            /*
                         * If the solar altitude is less than or equal to zero is night, so
                         * give up and go home, otherwise compute absolute optical air mass.
             */
            if (solel <= 0) {
                swave = 0;
            } else {
                rlpath = (Math.pow((solsin
                        + 0.15
                        * Math.pow(solel
                                + 3.885, -1.253)), -1));

                /* if(rlpath<1.0) rlpath=1.0; */
                path = (0.001 * ps1 * rlpath);
                solel = solel / radian;

                /**
                 * *******************************************
                 * transm ******************************************
                 */

                /*
                                 * transm calculates the solar transmission using the three-way
                                 * lookup table produced in gettbl
                 */
                ftabs = _transmft(rlpath, abstbl, ps1);
                ftscat = _transmft(rlpath, scatbl, ps1);
                fbscat = _transmfb(rlpath, bsctbl);

                /* Store values for use */
                tabs = ftabs;
                tscat = ftscat;
                bscat = fbscat;

                /* Set the path length for diffuse radiation equal to 1.7 */
                path = 1.7;

                /* Call transm again */
                ftabs = _transmft(path, abstbl, ps1);
                ftscat = _transmft(path, scatbl, ps1);
                fbscat = _transmfb(path, bsctbl);

                /* Again, store values for use */
                tabsd = ftabs;
                tscatd = ftscat;
                bscatd = fbscat;

                /*
                                 * sheat >>> sunlight amount on horizontal plane outside atmosphere
                                 * xser  >>> Diffuse shortware radiation at the ground.
                                 * hi    >>> Direct shortwave radiation reaching the ground
                                 * swave >>> Diffuse+direct = "Global"
                 */
                sheat = skonst * sinsolelslope / s;
                xser = bscatd * albdoe * (1 - tscatd) * tabsd
                        * Math.sin(solel);
                hi = (sheat * tabs * tscat)
                        + (skonst / s * tabs * (1 - tscat) * (1 - bscat))
                        * Math.sin(solel);
                swave = (hi * (1 - albdoe)) / (1 - xser);

                if (cloud_flag == 'T') {

                    /* Impose cloud fraction; reduce swave accordingly. */
                    swave = (swave * (1 - (0.7 * (0.01 * cld_fract))));
                }
            }

            /**
             * **********************************
             * end of input, back to netrad *********************************
             */

            /*
                         * Calcluate the net radiation for the appropriate ground conditions.
                         * the decision is made on the fraction of vegetation.
             */

 /*
                         * Code added 28/11/90 to solve the partial routine issue.
                         * Initialise and calculate the effective emissivity of the air and
                         * longwave down using a weighted average of surface and air temperatures.
             */
            if (init == 1) {
                aepsi = (0.7 + 0.17 * log10(omega));

                if (cloud_flag == 'T') {
                    aepsi = (aepsi + (1 - aepsi) * 0.8 * (0.01 * cld_fract));
                }

                bareradiotemp = tscren - 2;
                vegnradiotemp = tscren - 2;
                init = 2;
            }

            /* Vegetation */
            if (frveg == 1) {

                /**
                 * *************
                 * lwdown ************
                 */
                lwdn = aepsi * sigma
                        * Math.pow((t_fine[3] - tdif_s - 1.5), 4);

                /**
                 * *************
                 * end lwdown ************
                 */
                /**
                 * *********************
                 * vegrad ********************
                 */
                if (time == 0.0) {

                    /* initialise */
                    taf = otemp;
                    tf = otemp;
                    tg = otemp;
                    t1 = atemp;
                    ta = atemp;

                    /* qa=qd[1] */
                    qaf = qd[1];
                }

                /* Calculate incident solar flux at top of the canopy (sol) */
                sol = swave / (1 - albdoe);
                rsg = sol * (1 - sigf) * (1 - albg)
                        / (1 - sigf * albg * albf);
                rsf = sol * (1 - albf) * sigf
                        * (1 + albg * (1 - sigf) / (1 - sigf * albf * albg));
                rlg = (1 - sigf) * epsi * (lwdn - sigma * Math.pow(tg, 4))
                        / (1 - sigf * (1 - epsf) * (1 - epsi))
                        - epsi * epsf * sigf * sigma
                        * (Math.pow(tg, 4) - Math.pow(tf, 4))
                        / (1 - sigf * (1 - epsf) * (1 - epsi));
                rlf
                        = sigf
                        * (epsf * (lwdn - sigma * Math.pow(tf, 4)) + epsf * epsi * sigma * (Math.pow(tg, 4) - Math.pow(tf, 4)) / (1 - sigf * (1 - epsf) * (1 - epsi)))
                        + sigf * (1 - sigf) * (1 - epsi) * (epsf)
                        * (lwdn - sigma * Math.pow(tf, 4))
                        / (1 - sigf * (1 - epsi) * (1 - epsf));

                /*
                                 * Note that the radiometric temperature at the surface is denoted
                                 * by tzero and is the equivalent to otemp in the bare soil mode.
                 */
                tzero5 = (1 - sigf) * epsi * Math.pow(tg, 4)
                        + epsf * sigf * Math.pow(tf, 4);
                tzero5 = Math.pow(tzero5, 0.25);
                tzero4 = lwdn - rlg - rlf;
                tzero4 = tzero4 / (1 * sigma);
                vegnradiotemp = Math.pow(tzero4, 0.25);
                /* would be tzero */
                rnetg = rlg + rsg;
                rnetf = rlf + rsf;
                vegnnetradn = rnetg + rnetf;
                /* Would be rnetv */

 /* Set temperature and humidity variables */
                if (rnetf <= 0.0) {

                    /*
                                         * Set value since ta not set to t1 (vel.)
                                         * set in vel.
                     */
                    taf = ta;
                    tf = otemp;
                    tg = otemp;

                    /* ta=atemp; */
 /* qa=qd[1]; */
                    qaf = qd[1];
                }

                vegnshortwave = (rsg + rsf);
                /* Would be swavev */

                /**
                 * ***************
                 * end vegrad **************
                 */
                rnet = vegnnetradn;
                swave = vegnshortwave;
            } else if (frveg == 0) /* Bare soil */ {

                /* lwdown again */
                lwdn = aepsi * sigma
                        * Math.pow((t_fine[3] - tdif_s - 1.5), 4);

                /* uplong */
                lwup = epsi * sigma * Math.pow(bareradiotemp, 4);
                barenetradn = lwdn * epsi + swave - lwup;
                rnet = barenetradn;
            } else if ((frveg > 0.0) && (frveg < 1)) /* mixed */ {

                /* lwdown */
                lwdn = aepsi * sigma
                        * Math.pow((t_fine[3] - tdif_s - 1.5), 4);

                /* uplong */
                lwup = epsi * sigma * Math.pow(bareradiotemp, 4);
                barenetradn = lwdn * epsi + swave - lwup;

                /**
                 * **********
                 * vegrad **********
                 */
                if (time == 0.0) {

                    /* initialise */
                    taf = otemp;
                    tf = otemp;
                    tg = otemp;
                    t1 = atemp;
                    ta = atemp;

                    /* qa=qd[1] */
                    qaf = qd[1];
                }

                /* Calculate incident solar flux at top of the canopy (sol) */
                sol = swave / (1 - albdoe);
                rsg = sol * (1 - sigf) * (1 - albg)
                        / (1 - sigf * albg * albf);
                rsf = sol * (1 - albf) * sigf
                        * (1 + albg * (1 - sigf) / (1 - sigf * albf * albg));
                rlg = (1 - sigf) * epsi * (lwdn - sigma * Math.pow(tg, 4))
                        / (1 - sigf * (1 - epsf) * (1 - epsi))
                        - epsi * epsf * sigf * sigma
                        * (Math.pow(tg, 4) - Math.pow(tf, 4))
                        / (1 - sigf * (1 - epsf) * (1 - epsi));
                rlf
                        = sigf
                        * (epsf * (lwdn - sigma * Math.pow(tf, 4)) + epsf * epsi * sigma * (Math.pow(tg, 4) - Math.pow(tf, 4)) / (1 - sigf * (1 - epsf) * (1 - epsi)))
                        + sigf * (1 - sigf) * (1 - epsi) * (epsf)
                        * (lwdn - sigma * Math.pow(tf, 4))
                        / (1 - sigf * (1 - epsi) * (1 - epsf));

                /*
                                 * Note that the radiometric temperature at the surface is denoted
                                 * by tzero and is the equivalent to otemp in the bare soil mode.
                 */
                tzero5 = (1 - sigf) * epsi * Math.pow(tg, 4)
                        + epsf * sigf * Math.pow(tf, 4);
                tzero5 = Math.pow(tzero5, 0.25);
                tzero4 = lwdn - rlg - rlf;
                tzero4 = tzero4 / (1 * sigma);
                vegnradiotemp = Math.pow(tzero4, 0.25);
                /* would be tzero */
                rnetg = rlg + rsg;
                rnetf = rlf + rsf;
                vegnnetradn = rnetg + rnetf;
                /* Would be rnetv */

 /* Set temperature and humidity variables */
                if (rnetf <= 0.0) {

                    /*
                                         * Set value since ta not set to t1 (vel.)
                                         * set in vel.
                     */
                    taf = ta;
                    tf = otemp;
                    tg = otemp;

                    /* ta=atemp; */
 /* qa=qd[1]; */
                    qaf = qd[1];
                }

                vegnshortwave = (rsg + rsf);
                /* Would be swavev */

                /**
                 * ***************
                 * end vegrad **************
                 */
                mixednetradn = (vegnnetradn * frveg)
                        + (1 - frveg) * barenetradn;
                rnet = mixednetradn;
                swave = vegnshortwave * frveg + swave * (1 - frveg);
            } else {

                // System.out.println("Crapped out, frveg invalid = " + frveg + "\n");
                throw new Exception("Crapped out, frveg invalid = %f\n"
                        + frveg);
            }

            //frveg = frveg * 100;
            /**
             * ******************************************
             * end netrad *****************************************
             */

            /*
                         * Resistance values in the transition and surface layers.
                         * entry to nighttime formulations (bri and mom) through
                         * this subroutine.
             */
            /**
             * ******************************************
             * vel *****************************************
             */

            /*
                         * Vel computes the Monin Obukhov length, the friction velocity, and
                         * the integral of heat diffusivity. *?
                         * 
                         * Code altered 5th may 1992... Code transfrerred to bri.for
                         * alterations 16th july to account for different roughness lengths
                         * associated with partial vegetation calculations.
             */
 /* zo roughness height, za top of surface layer. (50m.) */
 /* zten - height at 10m, reflev -"Screen or anemometer height." */
            if (vel_pass == 0) {
                reflev = 2.0;
                zten = 10.0;
                ks = 0.0249951;
                kw = 0.0297681;
                cmh = 1;
                cmw = 1;
                vel_pass = 1;
            }

            /*
                         * The model assumes neutral conditions at the start of the run where
                         * heat=0. Therefore calc surface wind profile and resistances for the
                         * surface layer on the basis Math.log wind profile law.
             */
            if ((stabcriteria == 0)
                    || ((stabcriteria == 1) && (heat <= 0.0))) {

                /* Neutral profile - rare */
 /* Allow for negative heat flux */
                ustar = _you_star(awind, za, zo, 0.0);
                /* friction velocity */
                uscrn = _wind(ustar, reflev, zo, 0.0);
                /* wind speeds */
                uten = _wind(ustar, zten, zo, 0.0);
                rzazo = _r_ohms(ustar, za, zo, 0.0);
                /* Resistances */
                rzascr = _r_ohms(ustar, za, reflev, 0.0);

                if (dual_regime == 1) {
                    u_patch = _wind(ustar, obst_hgt, zo, 0.0);
                    ustar_patch = _you_star(u_patch, obst_hgt, zo_patch, 0.0);
                    rza_obst = _r_ohms(ustar, za, obst_hgt, 0.0);
                    robst_patch = _r_ohms(ustar_patch, obst_hgt, zo_patch,
                            0.0);
                }

                /*
                                 * Potential, actual temp, specific humidity at 2m. Pass
                                 * to vegetation component.
                 */
                ptmp20 = aptemp + (heat * rzascr / (DENS * CP));
                ta = ptmp20 - tdif_s;
                qa = ahum + (evap * rzascr / (DENS * le));
            } else if ((stabcriteria == 1) && (heat > 0)) {

                /*
                                 * Unstable case... use the non-dimensional functions of z/l in
                                 * wind profiles for momentum and heat (fm, ft) to determine the
                                 * resistance term over the surface layer.
                 */
                if (evap < 0.000001) {
                    evap = 0.000001;
                }

                bowen = heat / evap;

                if (Math.abs(bowen) > 10) {
                    bowen = 10 * bowen / Math.abs(bowen);
                }

                mol = (-1) * ((Math.pow(ustar, 3)) * aptemp * CP * DENS)
                        / (KARMAN * GRAV * heat
                        * (1 + 0.07 / (Math.abs(bowen))));

                /* Dimensionless wind shear */
                sa = _stab(za, mol);
                srfczo = _stab(zo, mol);
                so10m = _stab(reflev, mol);
                sten = _stab(zten, mol);

                /* Stability correction for momenttum. Benoit solution. */
                fm = _fstabm(srfczo, sa);
                ften = _fstabm(srfczo, sten);
                ustar1 = _you_star(awind, za, zo, fm);
                /* Friction velocity */
                ustar = (ustar + ustar1) / 2;
                /* Smooth */

 /* Dimensionless temperature gradient using the integrated form. */
                chio = _stabh(zo, mol);
                chia = _stabh(za, mol);
                chi20 = _stabh(reflev, mol);

                /* Stability correction for heat. Benoit solution */
                t_ft = _fstabh(chio, chia);
                ft20 = _fstabh(chi20, chia);
                uten = _wind(ustar, zten, zo, ften);
                /* Surface winds at */
                uscrn = _wind(ustar, reflev, zo,
                        0.0);
                /* screen level and 10m */
                rzascr = _r_ohms(ustar, za, reflev, ft20);
                /* Resistances in */
                rzazo = _r_ohms(ustar, za, zo, t_ft);
                /* surface layer */

                if (dual_regime == 1) {
                    sobst_hgt = _stab(obst_hgt, mol);
                    fobst_hgt = _fstabm(srfczo, sobst_hgt);
                    u_patch = _wind(ustar, obst_hgt, zo, 0.0);
                    ustar_patch = _you_star(u_patch, obst_hgt, zo_patch,
                            fobst_hgt);
                    ch_obst_hgt = _stabh(obst_hgt, mol);
                    ft_obst_hgt = _fstabh(ch_obst_hgt, chia);
                    rza_obst_hgt = _r_ohms(ustar, za, obst_hgt, ft_obst_hgt);
                    robst_patch = _r_ohms(ustar_patch, obst_hgt, zo_patch,
                            0.0);
                }

                /* As in neutral case, calculate TA and QA */
                ptmp20 = aptemp + (heat * rzascr / (DENS * CP));
                ta = ptmp20 - tdif_s;
                qa = ahum + (evap * rzascr / (DENS * le));
            } else {

                /* Call the nighttime routine. Stable case. */
                /**
                 * ********************************
                 * bri *******************************
                 */

                /*
                                 * Calculates the M-O-L when bri is (+) and less than 0.2
                                 * using the Blackadar model
                 */
 /* Code altered 5th may 1992...transfew from vel.for */

 /*
                                 * zan is the night time surface depth and is input in start.
                                 * zan is simply za at night.
                 */
                phih = 0;
                x1 = 0;

                if ((ycount <= 0) && (ifirst == 0)) {

                    /* mol=10e5; */
                    tstar = 0;
                    t[1] = aptemp;
                    t1 = tscren;
                    ifirst = 1;
                } else if ((ycount >= 1) && (ifirst == 1)) {

                    /* Start the surface temperature at some nominal value */
                    otemp = atemp;
                    t1 = atemp - 1;
                    tstar = 0;

                    for (i = 1; i <= 46; i++) {
                        t[i] = t_fine[i];
                        u[i] = u_fine[i];
                        v[i] = v_fine[i];
                        qn[i] = q_fine[i];
                    }

                    ifirst = 2;
                }

                /* Surface potential temperature */
                pot_s = t1 + tdif_s;

                /*
                                 * Calc the windspeed at the first level and the critical Richardson
                                 * number
                 */
                wg1 = Math.sqrt(ug[1] * ug[1] + vg[1] * vg[1]);
                cr1 = (Math.exp(-0.2129 * wg1) * 0.5542) + 0.2;

                do {
                    tdif = Math.abs(t[1] - pot_s);

                    /*
                                         * Calc the bulk Richardson number using the Blackadar scheme along
                                         * with the parameterisation for d0/dt
                     */
                    atemp = t[1] - tdif_50;
                    radcor = A * (otemp - t1) - rad;

                    /* tdif=Math.abs(t[1]-t1); */
                    tsurf = t1 - tstar * Math.log(Z1 / zo);
                    bulk = ((t[1] - pot_s) * GRAV * za)
                            / (otemp * awind * awind);

                    /*
                                         * Now use this to determine the stability criteria and execute the
                                         * appropriate physics
                     */
                    if (tdif < 0.05) {

                        /* Soln. sequence neutral */
                        tstar = (t[1] - pot_s) / Math.log(za / Z1);
                        ustar1 = KARMAN * awind / Math.log(za / zo);
                        ustar = (ustar1 + ustar) / 2;
                        heat = (-1) * KARMAN * DENS * CP * ustar * tstar;
                        uten = ustar / KARMAN * (Math.log(zten / zo));
                    } else if (bulk < 0) {

                        /* Soln. sequence unstable */

 /*
                                 * On entering this routine the mol will be positive wich will
                                 * result in a domain error. Hence restrain mol first time through.
                         */

 /*
                                 * if(mol>=0)
                                 * {
                                 * mol=-100000;
                                 * }
                                 * x=Math.pow((1-16*(za/mol)),-0.25);
                                 * phih=2*Math.log((1+x*x)/2);
                                 * phim=(phih+PI)/2+2*Math.log((1+1)/2)-2*Math.atan(x);
                                 * mol1=bulk/za*Math.pow(Math.log(za/z1)-phim,2)/(Math.log(za/z1)-phih);
                                 * mol=1/mol1;
                                 * mol=1e5;
                         */
                        ustar = KARMAN * awind / (Math.log(za / zo));
                        tstar = (t[1] - pot_s) / (Math.log(za / Z1));
                        heat = (-1) * KARMAN * DENS * CP * ustar * tstar;
                        uten = ustar / KARMAN * (Math.log(zten / zo));
                    } else if ((bulk > 0) && (bulk < cr1)) {

                        /* Solution sequence stable turbulent */
                        mol1 = (1 / za) * Math.log(za / zo)
                                * (bulk / (5 * (cr1 - bulk)));
                        mol = 1 / mol1;
                        ztenoverl = zten / mol;
                        psimzten = ANEW * ztenoverl
                                + BNEW * (ztenoverl - CNEW / DNEW)
                                * Math.exp(-DNEW * ztenoverl)
                                + BNEW
                                * CNEW
                                / DNEW;

                        /* A bunch of commented out stuff here */

 /*
                                 * cc	PSIHzten  =  (1.0 + 0.6667 * ANEW * ZtenOVERL)**1.5 + BNEW *
                                 * cc     &	   ( ZtenOVERL - CNEW / DNEW) * Math.exp ( - DNEW * ZtenOVERL) +
                                 * cc     &	   BNEW * CNEW / DNEW - 1.0
                                 * 
                                 * cc       ZrefOVERL = reflev / mol
                                 * 
                                 * cc	PSIMref = ANEW * ZrefOVERL + BNEW * (ZrefOVERL - CNEW / DNEW) *
                                 * cc	&	Math.exp ( - DNEW * ZrefOVERL) + BNEW * CNEW / DNEW
                                 * 
                                 * cc	   PSIHref = (1.0 + 0.6667 * ANEW * ZrefOVERL)**1.5 + BNEW *
                                 * cc	&	( ZrefOVERL - CNEW / DNEW) * Math.exp ( - DNEW * ZrefOVERL) +
                                 * cc	&	   BNEW * CNEW / DNEW - 1.0
                                 * 
                                 * cc	USCRN  =  USTAR / KARMAN * ( ALOG ( REFLEV / ZO ) + psimref )
                         */
 /* Compute the static stability correction */
                        zoverl = za / mol;
                        psimnew = ANEW * zoverl
                                + BNEW * (zoverl - CNEW / DNEW)
                                * Math.exp(-DNEW * zoverl)
                                + BNEW * CNEW
                                / DNEW;
                        psihnew = Math.pow((1 + 0.6667 * ANEW * zoverl), 1.5)
                                + BNEW * (zoverl - CNEW / DNEW)
                                * Math.exp(-DNEW * zoverl)
                                + BNEW * CNEW
                                / DNEW - 1;
                        tstar = (t[1] - pot_s)
                                / (Math.log(za / Z1) + psihnew);
                        ustar1 = KARMAN * awind
                                / (Math.log(za / zo) + psimnew);
                        ustar = (ustar1 + ustar) / 2;
                        heat = (-1) * KARMAN * DENS * CP * ustar * tstar;
                        uten = ustar / KARMAN
                                * (Math.log(zten / zo) + psimzten);
                    } else {

                        /* Soln sequence stable non-turbulent */
                        heat = -0.001;
                        ustar = ustar / 2;
                        uten = uten / 2;
                    }

                    if (ustar < 0.01) {
                        ustar = 0.01;
                    }

                    t2 = t1
                            + DT
                            * (radcor + advgt - B * heat / (CP * DENS * Z1));
                    t1 = (t1 + t2) / 2;
                    t1 = t1 - 0.017;

                    /* Mitre the night-time loop; cycle through twice (120s) */
                    x1 = x1 + 1;
                    xmod = x1 % 2;		  // was cast to (int)

                    /*
                                         * Fucking fucking stupid goto here originally. Stuffed everything
                                         * from where the label was up to here into a "do" loop.
                     */
                } while (xmod != 0);

                /*
                                 * Call mom to calculate the night-time vertical profiles of
                                 * temperature and winds
                 */
                if (ycount >= 1) {

                    /**
                     * *********************
                     * mom *********************
                     */

                    /*
                                         * mom calculates the momentum and thermodynamic eqns in the
                                         * lowest 500m of the atmosphere and produces profiles of the u and v
                                         * components, humidity and temperature at 50m intervals to couple
                                         * surface and mixing layers.
                     */
                    if (mom_pass == 0) {
                        krad = 0.75;

                        /*
                                 * The following were initialised inside a conditional, which
                                 * is meaningless since fortran evaluates them at compile time,
                                 * irrespective of the conditions.
                         */
                        dz = 50;
                        sb = 784;
                        imax = 46;
                        mom_pass = 1;
                    }

                    okmax = 1;
                    x1 = 0;

                    if (monce == 0) {

                        /*
                                 * Perform initialisation.
                                 * uif, uif2 are correction factors to account for 1st layer not
                                 * being centred when the night slf depth is below 50m (za.)
                                 * This is not necessary at present as za=50m.
                         */
                        uif_mom = (2 * dz / za) - 1;
                        uif2 = dz / (2 * dz - za);

                        /* Set up constants for the Richardson number calculation */
                        rc = (GRAV * dz) / otemp;
                        imax1 = imax - 1;

                        /*
                                 * Set values of critical Richardson # and turbulent
                                 * diffusivities to a nominal value and smooth the temperature
                                 * profile
                         */
                        for (i = 1; i <= imax1; i++) {
                            cr[i] = 0.001;
                            ok[i] = 0.01;
                            t[i] = (t[i] + t[i + 1]) / 2;
                        }

                        /* Make the diffusion coefficient (bx1) for radiation unitless */
                        bx1 = krad * DT / (dz * dz);
                        monce = 2;
                    }

                    /*
                                         * Calculate bulk parameters... windspeed, geostrophic windspeed
                                         * and critical Richardson number at all levels
                     */
                    do /* This was where the label for the goto was */ {
                        for (i = 1; i <= imax; i++) {
                            wind[i] = Math.sqrt(u[i] * u[i] + v[i] * v[i]);
                            wgeos[i] = Math.sqrt(ug[i] * ug[i]
                                    + vg[i] * vg[i]);
                            cr[i] = (Math.exp(-0.2129 * wgeos[i]) * 0.5542)
                                    + 0.2;
                        }

                        /*
                                 * Take derivatives of u and v and their absolute value. Set a
                                 * minimum value for the derivatives. Calculate the local Richardson
                                 * number using the expression for si and out constant (g/ ).
                                 * Next, find the values for the eddy exchange coefficients (kh & km);
                                 * Note... assumed to be the same in the stable nocturnal boundary
                                 * layer (sb=1). Smooth and set the non-dimensional forms for
                                 * momentum (b) and heat (bx.)
                         */
                        for (i = 2; i <= imax; i++) {
                            du = u[i] - u[i - 1];
                            dv = v[i] - v[i - 1];
                            abdu = Math.abs(du);
                            abdv = Math.abs(dv);

                            if (abdu < 0.001) {
                                du = 0.001;
                            }

                            if (abdv < 0.001) {
                                dv = 0.001;
                            }

                            ri[i] = rc * (t[i] - t[i - 1])
                                    / (du * du + dv * dv);
                            ok[i] = sb * ((cr[i] - ri[i]) / cr[i])
                                    * (Math.sqrt(du * du + dv * dv) / dz);
                            ok[i] = ok[i] + 0.05 * ok[i - 1];

                            if (ri[i] >= cr[i]) {
                                ok[i] = 0;
                            }

                            b_mom[i] = (ok[i] * DT) / (dz * dz);
                            bx[i] = ((ok[i] + krad) * DT) / (dz * dz);

                            if (b_mom[i] > 0.25) {
                                b_mom[i] = 0.25;
                            }

                            if (bx[i] > 0.25) {
                                bx[i] = 0.25;
                            }

                            if (ok[i] > okmax) {
                                okmax = ok[i];
                            }
                        }

                        /*
                                 * Calculate the vertical profiles of temperature and wind from 50
                                 * to 500m by integrating the u and v momentum equations and the
                                 * thermodynamic equation.
                         */
                        for (i = 2; i <= imax1; i++) {
                            if (i > 2) {
                                uif2 = 1;
                            }

                            duw[i] = b_mom[i + 1] * (u[i + 1] - u[i])
                                    - b_mom[i] * (u[i] - u[i - 1]) * uif2;
                            dvw[i] = b_mom[i + 1] * (v[i + 1] - v[i])
                                    - b_mom[i] * (v[i] - v[i - 1]) * uif2;
                            dtw[i] = bx[i + 1] * (t[i + 1] - t[i])
                                    - bx[i] * (t[i] - t[i - 1]) * uif2;
                            qn[i] = qn[i]
                                    + (b_mom[i + 1] * (qn[i + 1] - qn[i])
                                    - b_mom[i] * (qn[i] - qn[i - 1]));
                            u[i] = u[i] + cf * DT * (v[i] - vg[i]) + duw[i];
                            v[i] = v[i] - cf * DT * (u[i] - ug[i]) + dvw[i];
                            t[i] = t[i] + dtw[i] - (rad - advgt) * DT;
                        }

                        /*
                                 * Integrate at top and bottom boundaries where no turbulent
                                 * exchange exists and at the surface boundary condition.
                         */
                        t[1] = t[1] + heat * DT / (CP * DENS * za)
                                + bx[2] * (t[2] - t[1]) - (rad - advgt) * DT;
                        t[1] = t[1] - bx1 * (t[1] - pot_s);
                        u[1] = u[1] + cf * DT * (v[1] - vg[1])
                                + b_mom[2] * (u[2] - u[1])
                                - DT / dz * ustar * ustar * u[1] / wind[1]
                                * uif;
                        v[1] = v[1] - cf * DT * (u[1] - ug[1])
                                + b_mom[2] * (v[2] - v[1])
                                - DT / dz * ustar * ustar * v[1] / wind[1]
                                * uif;
                        qn[1] = qn[1]
                                + (b_mom[2] * (qn[2] - qn[1])
                                + (evap * DT / (DENS * le * dz)));
                        u[imax] = u[imax] + cf * DT * (v[imax] - vg[imax])
                                - b_mom[imax] * (u[imax] - u[imax1]);
                        v[imax] = v[imax] - cf * DT * (u[imax] - ug[imax])
                                - b_mom[imax] * (v[imax] - v[imax1]);
                        t[imax] = t[imax] - bx[imax] * (t[imax] - t[imax1])
                                - (rad - advgt) * DT;

                        /* Here we set qd=qn to be used in flux */
                        qd[1] = qn[1];

                        /* Increment the time control and cycle through twice. */
 /* And with another fucking goto, no less. But not any more. */
                        x1 = x1 + 1;
                        xmod = x1 % 2;			  // was cast to (int)
                    } while (xmod != 0);

                    awind = (Math.sqrt(u[1] * u[1] + v[1] * v[1]));

                    /* Fine mesh */
                    for (i = 1; i <= imax; i++) {
                        u_fine[i] = u[i];
                        v_fine[i] = v[i];
                        t_fine[i] = t[i];
                        q_fine[i] = qn[i];
                    }

                    /**
                     * ****************************
                     * end of mom ***************************
                     */
                }

                /**
                 * *******************************
                 * end of bri ******************************
                 */
                rekust = 1 / (KARMAN * ustar);
                rzazo = rekust * (Math.log(za / zo) + psihnew);

                /*
                                 * Note that stability correction function corresponds to
                                 * old method described in Panofsky and Dutton. New function
                                 * emplyed in bri,for. Differecnes between two are slight during
                                 * the short period following radiation sunset when the above
                                 * functions are employed.
                 */
                ta = t1;
                qa = ahum;
            }

            /* Set ustar equal to some non-zero value if small */
            if (ustar < 0.01) {
                ustar = 0.01;
            }

            /*
                         * Calc the diffusivities and resistances to heat and water for
                         * the transition layer. ks and kw are the molecular conductivities
                         * for heat and water. cmh and cmw are the scaling factors. See manual
             */
            if (ionce == 0) {
                kmm = cmh * ks / (DENS * CP);
                kx = cmw * kw / (DENS * CP);
            }

            rtrans = _restrn(ustar, zo, kmm);
            /* resistance transition layer */
            rtranw = _restrn(ustar, zo, kx);
            rtrano3 = _restrn(ustar, zo, kx / 1.32);

            /* Finally, compute the total series resistances over both layers */
            sum = rzazo + rtrans;
            sumw = rzazo + rtranw;
            sumo3 = rzazo + rtrano3;

            if (dual_regime == 1) {
                rtrans_patch = _restrn(ustar_patch, zo, kmm);
                rtransw_patch = _restrn(ustar_patch, zo, kx);
                sum = rza_obst_hgt + robst_patch + rtrans_patch;
                sumw = rza_obst_hgt + robst_patch + rtransw_patch;
            }

            /* If vegetation is included, call the vegetation component */
            if ((swave > 0) && (rnet > 0) && (frveg > 0)) {

                /**
                 * *************************************
                 * vegvel ************************************
                 */
                if (vegvel_pass == 0) {
                    init_vel = 1;
                    vegvel_pass = 1;
                }

                pes = (xlai / 2.0) + 1;
                rhocp = CP * ps1 * 100 / (r * taf);
                sigalf = 1 - 0.5 / (0.5 + xlai) * Math.exp(-xlai * xlai / 8);

                /*
                                 * cdl is the drag coefficient for the lef.
                                 * Kel's formulation.
                 */
                cdl = 0.08;
                sdl = xlai / vegheight;
                uaf = ustar * Math.sqrt(pes / (cdl * sdl));
                chf = 0.011 * Math.sqrt(uaf / width) * (xlai / pes);

                if (chf < 0.001) {
                    chf = 0.001;
                }

                raf = 1 / chf;

                /*
                                 * c	CHG = ( 1 - SIGALF ) * USTAR**2 / UAF
                                 * c        Kell's chg
                 */
                chg = KARMAN * KARMAN * uaf
                        / Math.pow((Math.log(vegheight * 100 / 0.7)), 2);

                /*
                                 * Add surface layer resistance to internal air resistance and
                                 * recompute conductance
                 */
                rhg = 1 / chg + rtranw;
                chg = 1 / rhg;

                /*
                                 * 
                                 * C     FOLLOWING IS KEL'S FORMULATION FROM GOUDRIAAN
                                 * C     TO GET CHA
                                 * 
                                 * c      ALEN = 2 * ((3 * WIDTH**2) / (4 * PI * XLAI/(VEGHEIGHT/2)))
                                 * c     /	       ** (0.3333)
                                 * c     AKCAN = ALEN * 0.5 * UAF
                                 * c      RHA = VEGHEIGHT / (2 * AKCAN)
                                 * c      CHA = 1 / RHA
                                 * c       Above is outdated formulae from Kell. Use old method
                 */
                cha = ustar * ustar / (uten - uaf);

                if (cha <= 0.001) {
                    cha = 0.001;

                    /* Convert q to eltl and ea */
                }

                qstf = Math.pow(10, (6.1989 - (2353 / tf)));
                sgma = rhocp / 0.666;
                rmratiodif = qstf - qaf;

                if (rmratiodif <= 0) {
                    rmratiodif = 0.0001;
                }

                vfl = rmratiodif * ps1 / 0.622;

                /*
                                 * Call Deardorff (D) or Carlson/Lynn formulation for RST.
                                 * Note rst is actually rf because it contains a cuticular resistance,
                                 * rcut.
                 */
                if (stmtype == 'D') {
                    //System.out.println("-- Deardorff");
                    /* Deardorff formulation */
 /* Stomatal resistance */
                    rs = rmin
                            * (800.0 / (1.0 + sol)
                            + Math.pow((1.2 * wilt / (0.9 * w2g + 0.1 * wgg)),
                                    2));

                    /* Total leaf/canopy resistance */
                    rst = rs * rcut / (rs + rcut) * pes / xlai;
                } else if (stmtype == 'L') {

                    //System.out.println("-- Carlson & Lynn");
                    /* Lynn and Carlson */
                    /**
                     * ********************
                     * PSLCAL ********************
                     */
                    if (pslcal_pass == 0) {
                        init_pslcal = 1;
                        pslcal_pass = 1;
                    }

                    /* Subroutine that finds leaf water potential analytically */

 /*
                                         * 17th june 1992- b1_p and b2_p are defined differently from original
                                         * Lynn and Calrson article. To allow rmin to be a true scaling coeff
                                         * for the stomatal resistance.
                     */

 /*
                                         * Explanation of conversion factors *
                                         * gma=(p*cp)/(0.622*le)=0.66mb c-1
                                         * leef=roe*le*v/(rl+raf)
                                         * substitute for le
                                         * leef=(roe*cp)/gma*v/(rl+raf)
                                         * e=q*(p-e)/0.622=q*p/0.622
                                         * p is set equal to a ps1 (mbs) below.
                                         * note sensible flux=roe*cp*delta T/rs
                                         * note conversion to latent flux is roe*cp*delta e/(rs*gma)
                                         * roe=(kgm***3)
                                         * cp=1005 k/kg-k
                                         * set thv=w2g
                     */
                    thv = w2g;
                    fs = _stomfs(sc, sol);

                    /*
                                         * Following is Choudury and Idso's adaptation for g
                                         * units is a conversion factor from s given in Choudury to
                                         * bar/(w/m**2)
                     */
                    units2 = 0.4e-10;

                    /*
                                         * veghieght (h) is the height of the plants in meters
                                         * zeff is the effective rooting depth taken as half of h.
                     */
                    zeff = 0.5 * vegheight;
                    rkw = _cond(rks, thv, thmax, cosbyb);
                    zg = (0.0013 / (zeff * rkw)) * units2;
                    ztot = zp + zg;

                    /*
                                         * Compute the total soil-root and root-xylem resistance.
                                         * frht is fraction of height below the intersection point
                                         * of the storage area with the xylem.
                                         * fpm represents the hieght of the tree above the
                                         * intersection point *?
                                         * 
                                         * Following is the adaptation of fw using the weight of the column
                                         * of water after Federer (1982).
                                         * This makes an adjustment of the hiehgt of the column.
                                         * The density of water of 1000kg/m**3
                                         * the acceleration of gravity is 9.8m/s**2
                                         * we need the average height of the plant or the weighted mean
                                         * height which is the displacement or the geometric height. The
                                         * displacement height is set at 0.67 times the height of the
                                         * vegetation.
                                         * units1 is a conversion factor from kg/s**2 to bars.
                     */
                    b1_p = b1 * rmin;
                    b2_p = b2 * rmin;
                    hbar = vegheight * 0.67;
                    units1 = 1e-5;
                    h = RHOW * grav * hbar * units1;
                    unscaled_raf = raf * xlai / pes;
                    /* Per leaf area */
                    psig = _psgcal(thmax, thv, psis, cosbyb);
                    rscrit = _stomc(b1, psice, rmin, fs, FT);

                    if (init_pslcal == 1) {
                        vfl2 = vfl;
                        init_pslcal = 2;
                    } else {
                        q_boundary = (qstf / rst + qaf / raf)
                                / (1 / raf + 1 / rst);
                        vfl2 = (qstf - q_boundary) * ps1 / 0.622;
                    }

                    //steady = 'J';
                    if (steady == 'Y') /* Capacitance solution */ {

                        /**
                         * *********************
                         * capac *********************
                         */
                        zadd = zg + zp * frzp;
                        fpm = 1 - frhgt;
                        zuh = zp * (1 - frzp);
                        //p.println(zadd + " " + fpm + " " + zuh + " " + zg);
                        /*
                                 * Below there are two solutions for psist. The first is the
                                 * Math.log, while the second is a regular integral based on volume
                                 * to a power. There are also two similar solutions for zst
                                 * and capacitance.
                         */

 /*
                                 * Note also that psist, psist and psix are a time step behind
                                 * because psist is not allowed to exceep psix of the previous time
                                 * step. (psix is a function of psis.)
                         */

 /*
                                 * There are also ways to set the volume of water in
                                 * storage. If xcap is set equal to zero, then volist storage is a
                                 * function of psig. If not, it is a function of xcap or
                                 * relative water content.
                         */
                        if (jcap == 1) {
                            ncap = (int) rccap;
                            idel = 0;
                            jdel = 0;
                            psix = psig - frhgt * h;
                            ixcap = (int) volrel;
                            volrel = 0.01 * volrel;
                            capini = volini / rkocap;
                            jcap = 2;

                            if (ixcap == 0) {
                                if (ncap == 1) {
                                    volist = volini
                                            * Math.exp(capini
                                                    * (psig - frhgt * h)
                                                    / volini);
                                    capaci = capini * (volist) / volini;
                                } else if (ncap != 1) {
                                    caprat = (rccap - 1) / rccap;
                                    volist
                                            = volini
                                            * Math.pow((1
                                                    + capini
                                                    * (psig - frhgt * h)
                                                    * caprat / volini), (1
                                                    / caprat));
                                    capaci = capini
                                            * Math.pow((volist / volini),
                                                    (1 / rccap));
                                }

                                /*
                                                 * It is possible that the volist will be calculated as
                                                 * being less than zero. This impossible. Therefore, we set:
                                 */
                                if (volist < 0) {
                                    volist = 0.00001 * volini;
                                    volrmv = 0.99999 * volini;
                                }
                            } else {
                                volist = volrel * volini;
                                voliso = volist;

                                if (ncap == 1) {
                                    capaci = capini * volrel;
                                } else if (ncap != 1) {
                                    caprat = (rccap - 1) / rccap;
                                    capaci = capini
                                            * Math.pow(volrel, (1 / rccap));
                                }
                            }
                        }

                        volrel = volist / volini;

                        if (ncap == 1) {
                            psist = (volini / capini) * Math.log(volrel);
                            capaci = capini * (volist) / volini;
                        } else {
                            psist = (volini / capini)
                                    * (Math.pow(volrel, caprat) - 1)
                                    * (rccap / (rccap - 1));

                            /* This statement is in place of volume constants */
                            if (psist >= 0) {
                                psist = 0;
                            }

                            /*
                                         * The following if statement causes psist to fall as volume
                                         * approaches zero (the abover equation reaches a limit.)
                                         * (itrap is defined below to prevent the loop from being
                                         * exercised on the first time step.)
                                         * The next statement prevents psist from becoming more
                                         * positive than psig
                             */
                            if ((itrap == 1) && (deltvst == 0)) {
                                psist = psix;
                            }

                            capaci = capini * Math.pow(volrel, (1 / rccap));
                        }

                        if (ncap == 1) {
                            zst = zstini / volrel;
                        } else if (ncap != 1) {
                            zst = zstini * Math.pow((1 / volrel), rzcap);
                        }

                        /*
                                 * The following statement prevents zst from becoming
                                 * very large
                         */
                        if (zst > 1e3) {
                            zst = 1e3;
                        }

                        //System.out.println("NONONONONONONONON");
                        /*
                                 * The following if statements prevent psist from exceeping
                                 * psix during a time interval
                         */
                        if ((deltvst < 0) && (psist < psix)) {
                            psist = psix;
                        } else if ((deltvst > 0) && (psist > psix)) {
                            psist = psix;
                        }

                        rstdiv = unscaled_raf
                                + (rcut * rscrit) / (rcut + rscrit);

                        /* We need to compare psig to psigc before proceeding */

 /*
                                 * We can determine this by solving for psig given various
                                 * parameters.
                         */
 /* Addit is an additional term that incorporates storage flux */
                        addit = 1 + zuh / zadd + zuh / zst;

                        /*
                                 * Define a critical water potential. As the flux from
                                 * storage approaches zero, it approaches that of the soil
                         */
                        psiwc = sgma_capac * vfl * zuh / rstdiv + psice
                                + beta * vfl2 + fpm * h;
                        psisup = psix;

                        if (psisup > psiwc) {
                            aroot = fs * FT * b1_p * (rcut + unscaled_raf)
                                    * (zst * (-1) - zadd);
                            broot
                                    = fs * FT * (rcut + unscaled_raf)
                                    * (rmin * (-zst - zadd) + b1_p * (zst * (psig - beta * vfl2 - h) + zadd * (psist - beta * vfl2 - fpm * h) - zadd * zst * sgma_capac * vfl * addit / (rcut + unscaled_raf)))
                                    + rcut * unscaled_raf * (-zst - zadd);
                            croot
                                    = fs * FT * (rcut + unscaled_raf)
                                    * (rmin * (zst * (psig - beta * vfl2 - h) + zadd * (psist - beta * vfl2 - fpm * h) - zadd * zst * sgma_capac * vfl * addit / (rcut + unscaled_raf)))
                                    + rcut * unscaled_raf
                                    * (zst * (psig - beta * vfl2 - h) + zadd * (psist - beta * vfl2 - fpm * h))
                                    - rcut * zadd * zst * sgma_capac * vfl
                                    * addit;
                        } else {
                            aroot = fs * FT * b2_p * (rcut + unscaled_raf)
                                    * (zst + zadd);
                            broot
                                    = fs * FT * (rcut + unscaled_raf)
                                    * ((rmin + b1_p * psice + b2_p * psice) * (-zst - zadd) - b2_p * (zst * (psig - beta * vfl2 - h) + zadd * (psist - beta * vfl2 - fpm * h) - zadd * zst * sgma_capac * vfl * addit / (rcut + unscaled_raf)))
                                    + rcut * unscaled_raf * (-zst - zadd);
                            croot = fs * FT * (rcut + unscaled_raf)
                                    * (rmin + b1_p * psice + b2_p * psice)
                                    * (zst * (psig - beta * vfl2 - h) + zadd * (psist - beta * vfl2 - fpm * h) - zadd * zst * sgma_capac * vfl * addit / (rcut + unscaled_raf))
                                    + rcut * unscaled_raf
                                    * (zst * (psig - beta * vfl2 - h) + zadd * (psist - beta * vfl2 - fpm * h))
                                    - rcut * zadd * zst * sgma_capac * vfl
                                    * addit;
                        }

                        actsqrt = Math.pow(broot, 2) - 4 * aroot * croot;
                        psie = (-broot - Math.sqrt(actsqrt)) / (2 * aroot);


                        /* p.println("aroot: " + aroot);
                                p.println("broot: " + broot);
                                p.println("croot: " + croot);
                                p.println("actsqrt: " + actsqrt);
                
                                p.println("*PSIE*= " + psie); 
                
                                p.println("frveg " + frveg);
                                p.println("fs " + fs);
                                p.println("FT " + FT);
                                p.println("b1_p " + b1_p);
                                p.println("rcut " + rcut);
                                p.println("unscaled_raf " + unscaled_raf);
                                p.println("zst " + zst);
                                p.println("zadd " + zadd);
                                p.println("rmin " + rmin);
                                p.println("psig " + psig);
                                p.println("beta " + beta);
                                p.println("vfl2 " + vfl2);
                                p.println("h " + h);
                                p.println("sgma_capac " + sgma_capac);
                                p.println("vfl " + vfl);
                                p.println("addit " + addit);
                                p.println("psist " + psist);
                                p.println("fpm " + fpm);
                                p.println("------------------------------"); */
                        //System.out.println("PSIE= " + psie);

                        /* stomrs again */

 /*
                                 * Calculate stomatal resistance coefficients are
                                 * initialised in stmcof.for
                         */
                        if ((tf > mintemp) && (tf < maxtemp)) {
                            if (psisup > psiwc) {
                                fpsie = 1 + b1 * psie;
                            } else {
                                fpsie = 1 + b1 * psice + b2 * (psice - psie);
                            }
                        } else {
                            rs = 5000;
                        }

                        rs = rmin * fs * fpsie * FT;

                        /* end stomrs again */
 /* Stress indices */
                        psim = psie + beta * vfl2;
                        pes = xlai / 2 + 1;
                        rlelf = sgma_capac * vfl
                                / (rs * rcut / (rs + rcut) + unscaled_raf);
                        psix = rlelf * zuh + psim + fpm * h;
                        fst = (-1) * (psist - psix) / zst;
                        deltvst = (-1) * delta * (psist - psix)
                                / (zst * le * RHOW);
                        fluxgd = (psig - psix - frhgt * h) / zadd;
                        volrmo = volrmv;
                        volist = volist + deltvst;
                        volrmv = volini - volist;

                        /*
                                 * The following statements prevent the volume in storage from
                                 * becoming negative
                         */
                        if ((volrmv > volini) && (deltvst <= 0)) {
                            if (idel == 0) {
                                deltvst = (volini - volrmo);
                                idel = 1;
                            } else {
                                deltvst = 0.0;
                                itrap = 1;
                            }

                            volrmv = 0.99999 * volini;
                            volist = 0.00001 * volini;
                        }

                        /**
                         * *********************
                         * end of capac *********************
                         */
                        //System.out.println("What am I doing here?");
                    } else {
                        rstdiv = unscaled_raf
                                + (rcut * rscrit) / (rcut + rscrit);
                        psiwc = sgma * vfl * ztot / rstdiv + psice
                                + beta * vfl2 + h;
                        psisup = psig;

                        if (psisup > psiwc) {
                            aroot = fs * FT * b1_p * (rcut + unscaled_raf)
                                    * (-1);
                            broot
                                    = fs * FT * (rcut + unscaled_raf)
                                    * (-rmin + b1_p * (psig - beta * vfl2 - h - ztot * sgma * vfl / (rcut + unscaled_raf)))
                                    - rcut * unscaled_raf;
                            croot
                                    = fs * FT * (rcut + unscaled_raf) * rmin
                                    * (psig - beta * vfl2 - h - ztot * sgma * vfl / (rcut + unscaled_raf))
                                    + rcut * unscaled_raf
                                    * (psig - beta * vfl2 - h)
                                    - rcut
                                    * ztot
                                    * sgma
                                    * vfl;
                        } else {
                            aroot = fs * FT * b2_p * (rcut + unscaled_raf);
                            broot
                                    = fs * FT * (rcut + unscaled_raf)
                                    * (-1 * (rmin + b1_p * psice + b2_p * psice) - b2_p * (psig - beta * vfl2 - h - ztot * sgma * vfl / (rcut + unscaled_raf)))
                                    - rcut * unscaled_raf;
                            croot = fs * FT * (rcut + unscaled_raf)
                                    * (rmin + b1_p * psice + b2_p * psice)
                                    * (psig - beta * vfl2 - h - ztot * sgma * vfl / (rcut + unscaled_raf))
                                    + rcut * unscaled_raf
                                    * (psig - beta * vfl2 - h)
                                    - rcut
                                    * ztot
                                    * sgma
                                    * vfl;
                        }

                        actsqrt = Math.pow(broot, 2) - 4 * aroot * croot;
                        psie = (-broot - Math.sqrt(actsqrt)) / (2 * aroot);
                        //p.println("PSIE2: " + psie);
                        psim = psie + beta * vfl2;

                        /**
                         * **********************
                         * stomrs *********************
                         */

                        /*
                                 * Calculate stomatal resistance coefficients are
                                 * initialised in stmcof.for
                         */
                        if ((tf > mintemp) && (tf < maxtemp)) {
                            if (psisup > psiwc) {
                                fpsie = 1 + b1 * psie;
                            } else {
                                fpsie = 1 + b1 * psice + b2 * (psice - psie);
                            }
                        } else {
                            rs = 5000;
                        }

                        rs = rmin * fs * fpsie * FT;

                        /**
                         * **********************
                         * end of stomrs *********************
                         */
                    }
                    /* end of steady block */

 /* Stress indices */
                    wpsi = psiwc - psisup;
                    rlpsi = psice - psie;

                    /**
                     * *****************
                     * end of pslcal *****************
                     */
                    rl = rs * rcut / (rs + rcut);
                    /* Leaf resistance */
                    rst = rl * pes / xlai;
                    /* Total leaf/canopy resistance */
                } else if (stmtype == 'B') {
                    System.out.println("-- Ball");

                    /**
                     * ***************
                     * ball ***************
                     */

                    /* 24th march 1995. rst is always unscaled in this routine */

 /*
                                         * A brand new routine to get fco2 using photo model.
                                         * 
                                         * par	>> par in mol units
                                         * rrtot	>> total resistance for co2 diffusion in mol units
                                         * rair	>> surface layer resistance
                                         * fco2	>> co2 flux (either scaled per leat or total
                                         * cii	>> internal leaf co2 in ppm
                                         * 
                                         * The following are from the two Farquhar papers and the
                                         * "Note on Carbon fixation modelling" contains the current
                                         * parameterizations for ccomp, ko, kc, wr. Most of the
                                         * parameters are now set according to Collatz, Ball, Berry.
                                         * 
                                         * akc >> mm constant for co2
                                         * ako >> mm constant for o2
                                         * oxy >> oxygen concentration in chloroplasts
                                         * vcmax >> maximum carboxylation rate
                                         * alpi >> quantum yield (photons per co2)
                                         * ccomp >> co2 compensation point
                                         * fco2 >> carbon exchange rate
                                         * resp >> dark respiration
                                         * cii >> internal co2 concentration
                                         * we >> rubisco limited photosynthesis
                                         * wr >> rubp regeneration limited photosynthesis
                                         * vc >> carbon exchange in absence of photorespiration
                                         * parq >> par when wr=we for given conditions
                                         * tau >> co2/o2 specificty ratio.
                                         * css >> leaf sfc co2 conductance
                                         * gs >> stomatal conductance
                     */

 /*
                                         * Carbon flux and resistances should be per unit leaf area
                                         * to get physiological cii
                     */
                    rst = rst * xlai / pes;
                    raf = raf * xlai / pes;
                    /* unscaled */
                    fco2 = fco2 * pes / xlai / frveg;
                    /* unscaled */

 /* par in mol units */
                    par = ((sol / 2) * 4.57) * 1e-6;
                    /* in mol units */
                    rair = rha + rzascr;
                    rrtot = 1.32 * raf + 1.66 * rst + rair;
                    rrtot = rrtot / 40.0;
                    /* in mol units */

 /* Set some parameters photosyn model (now read in start */
                    akco = 3.0e-4;
                    akoo = 0.300;
                    oxy = 0.209;
                    vcmaxo = 2.0e-4;
                    /* Collatz; 7.5e-5 data 1992 */
                    alpi = 12;
                    respo = 3.0e-6;
                    tauo = 2600;
                    ajpar = 380e-6;
                    ampar = 9.0;
                    /* conductance model */

 /*
                                         * Correct kinematic prooperties to temperature.
                                         * see Collatz, Ball, Berry, Farquhar
                     */
                    akc = akco * Math.pow(2.1, ((tf - 298) / 10));
                    tau = tauo * Math.pow(0.57, ((tf - 298) / 10));
                    ako = akoo * Math.pow(1.2, ((tf - 298) / 10));
                    vcmax = vcmaxo * Math.pow(2.4, ((tf - 298) / 10));
                    resp = respo * Math.pow(2.0, ((tf - 298) / 10));
                    ccomp = oxy / (2 * tau);
                    ajmax = ajpar
                            * (1.0 + 0.0409 * (tf - 303.0)
                            - 1.54e-3 * Math.pow(tf - 303.0, 2)
                            - 9.42e-5 * Math.pow(tf - 303.0, 3));

                    /*
                                         * Correct restp inhibition at high temperatures.
                                         * from Collatz, Ball, Berry
                     */
                    resp = resp * 1 / (1 + Math.exp((1.3 * (tf - 328))));

                    /*
                                         * Do iterations to get cii stable when conditions changing
                                         * rapidly
                     */
                    ball_repeat = 0;

                    do {

                        /* Guess at ci */
                        ci = co - fco2 * rrtot;
                        cii = co - fco2 * rrtot;
                        /* previous value */

                        if (ci < 0) {
                            ci = 220e-6;
                        }

                        /* Calculate new rubisco and rubp limited rates */
                        we = vcmax * ci / (ci + akc * (1 + (oxy / ako)));
                        aj = ajmax * par / (par + 2.1 * ajmax);
                        wr = aj * ci / (4.5 * ci + 10.5 * ccomp);
                        ws = vcmax / 2;

                        /* Determine which factor is limiting */
                        if ((we <= wr) && (we <= ws)) {
                            limit = 'R';
                            /* rubisco */
                        } else if ((wr < we) && (wr < ws)) {
                            limit = 'E';
                            /* E trans */
                        } else {
                            limit = 'S';
                            /* sink */
                        }

                        /* Set carboxylation rate to the limiting factor */
                        vc = _min(_min(we, wr), _min(wr, ws));
                        parq = we * (alpi * (1 + (2 * ccomp / ci)));

                        /* co2 flux */
                        fco2 = vc * (1 - ccomp / ci) - resp;

                        if (ccomp > ci) {
                            fco2 = (-1) * resp;
                        }

                        /* Ball Berry stomatal conductance model */
 /* gs=ampar*fco2*rh1/css+0.01; /* in css units */

 /*
                                 * Estimate of gs - recheck this: there may be a problem
                                 * using 1.32 for water vapor.
                         */
                        aaa = ci;
                        bbb = co / (raf * 1.32 + rair) - ampar * fco2;
                        ccc = (-1) * ampar * fco2 * qaf / (qstf * raf);

                        /*
                                 * Solution if no relative humidity dependence
                                 * following two statements
                         */
 /* bbb=co/(raf*1.32+rair); */
 /* ccc=(-1)*ampar*fco2; */
                        if ((bbb * bbb) < (4 * aaa * ccc)) {
                            rst = rcut;
                        } else {
                            gs = (-bbb + Math.sqrt(bbb * bbb - 4 * aaa * ccc))
                                    / (2 * aaa);
                            gs = gs / 40.0;
                            /* mks units */

                            if (gs != 0) {
                                rs = 1.0 / gs;
                            }

                            if (gs == 0) {
                                rs = 500;
                            }

                            rst = rs * rcut / (rs + rcut);
                            rrtot = 1.32 * raf + 1.66 * rst + rair;
                            rrtot = rrtot / 40.0;
                            /* mole units */
                            ci = co - rrtot * fco2;

                            /* Check if ci is converging */
                            if (Math.abs(1 - ci / cii) > 0.05) {
                                ball_repeat = 1;
                            }
                        }
                    } while (ball_repeat != 1);

                    if (fco2 < 0.0) {
                        rst = rcut;
                    }

                    /*
                                         * Scale the resistance and co2 flux. Not used in any
                                         * equations, just for output.
                     */
                    qbnd = (qstf / rst + qaf / raf)
                            * Math.pow((1 / raf + 1 / rst), -1);

                    /* humidity at leaf surface or leaf boundary layer */
                    rhl = qbnd / qstf;
                    /* rh at leaf surface */
                    ccan = (ci / (rst * 1.66 + raf * 1.32) + co / (rair))
                            * Math.pow(1 / (rst * 1.66 + raf * 1.32) + 1
                                    / rair, -1);
                    css = fco2 * raf * 1.66 + ccan;
                    /* co2 at leaft surface */
                    raf = raf * pes / xlai;
                    rst = rst * pes / xlai;
                    fco2 = fco2 * xlai / pes * frveg;

                    /**
                     * **********************
                     * end ball *********************
                     */
                    rl = rs * rcut / (rs + rcut);
                    rst = rl * pes / xlai;
                }

                qsta = Math.pow(10, (6.1989 - (2353 / ta)));
                xlefn = DENS * le * rmratiodif / (rst + raf);

                /* Average tf and xlef after first time step */
                if (init_vel == 1) {
                    tf = taf + (rnetf - xlefn) / (chf * DENS * CP);
                    xlef = xlefn;
                    init_vel = 2;
                } else {
                    xlef = (xlef + xlefn) / 2;
                    tfn = taf + (rnetf - xlefn) / (chf * DENS * CP);
                    tf = (tf + tfn) / 2;
                }

                /**
                 * *************************************
                 * end of vegvel ************************************
                 */
            }

            /**
             * *************************************
             * end of vel ************************************
             */

            /* Mixed layer */
            if ((heat > 0.00001) && (swave > 0) && (rnet > 0)) {

                /**
                 * ******************************************
                 * air *****************************************
                 */

                /*
                                 * air computes the daytime height of the mixing layer and
                                 * the potential temperature at height za
                 */
                ifirst = 1;

                /* signal daytime situation, ifirst=1 */
 /* Select the correct pot. temp. lapse rate */
                for (j = 2; j <= 9; j++) {
                    if (hgt > zls[j - 1]) {
                        gam = gm[j];
                    }
                }

                /*
                                 * chgt is calculated here once only, based initially on a
                                 * parameterisation of Tennekes, before passing onto the ELSE
                                 * statement. Now mitre the mixing layer loop to 24 sec to increase
                                 * accuracy, calc pot temp, temp at za and the mixing layer height
                                 * over 240 sec.
                 */
                if (chgt == 0) {
                    het = heat / (DENS * CP);
                    chgt = (0.35 * Math.sqrt(0.3))
                            * (Math.pow((grav / otemp), (1.0 / 3.0))
                            * Math.pow(hgt, (1.0 / 3.0))
                            * Math.pow(het, (1.0 / 3.0)));
                    cdelt = (gam * hgt * chgt - het) / hgt;
                    ctheta = (het / hgt) - rad + advgt;
                } else {
                    deltx = delta / 10;

                    for (i = 1; i <= 10; i++) {
                        aptemp = (ctheta * deltx) + aptemp;
                        atemp = aptemp - tdif_50;
                        delt = (cdelt * deltx) + delt;

                        if (delt < 0.01) {
                            delt = 0.01;
                        }

                        hgt = (chgt * deltx) + hgt;
                        het = heat / (DENS * CP);
                        dhet
                                = (-0.5) * het
                                / ((1.0
                                + (2.6 * Math.pow(het, (2.0 / 3.0))
                                / (Math.pow((grav * hgt / otemp), (1.0 / 3.0))
                                * delt))));
                        chgt = (-1) * dhet / delt;
                        ctheta = ((het - dhet) / hgt) - rad + advgt;
                        cdelt = ((gam * hgt * chgt) - (het) - (delt * chgt))
                                / hgt;
                    }
                }

                /* Sounding profile for course */
                td[1] = aptemp;
                zmix = za;

                for (i = 2; i <= ntrp; i++) {
                    zmix = zmix + deltaz;

                    if (zmix < hgt) {
                        td[1] = aptemp;
                    }
                }

                /* tdel is at the height just above the mixing layer */
 /* ttop is at the height of the mixing layer */
                tdel = aptemp + cdelt;
                ttop = aptemp;

                /* ycount advanced to set the mode indefinitely for the daytime. */
                ycount = ycount + 1;

                /**
                 * ******************************************
                 * end of air *****************************************
                 */
            }

            /* Eddy diffusivities in the mixed layer */
            if ((heat > 0.00001) && (swave > 0) && (rnet > 0)) {

                /**
                 * ******************************************
                 * daykm *****************************************
                 */
                thick_mx = _daykm(km, zi, zk, hgt, za, ustar, mol, ntrp,
                        deltaz);
            }

            /* Momentum equations - mixed layer */
            if ((heat > 0.00001) && (swave > 0) && (rnet > 0)) {

                /**
                 * ******************************************
                 * momday *****************************************
                 */

                /*
                                 * This daytime routine updates the daytime winds ud, vd and specific
                                 * humidity qd using the eddy diff's obtained in daykm, assuming
                                 * similarity between humidity and momentum transfer coefficients.
                                 * Note that the time step is 120 sec... this ensures computational
                                 * stability during periods when km values are large.
                                 * 
                                 * Initialise some variables
                 */
                if (momday_pass == 0) {
                    smfh = 60;
                    uif = 1;
                    momday_pass = 1;
                }

                /*
                                 * smf and thick is to deal with the eddy equations when the
                                 * mixing layer is below 60m
                 */
                smf = thick_mx / smfh;

                if (smf > 1) {
                    smf = 1;
                }

                awind = Math.sqrt(vd[1] * vd[1] + ud[1] * ud[1]);
                xxx = 0;

                /*
                                 * deltzi is deltaz for the average of the two lowest layers for
                                 * finite differencing of bottom boundary in diffusion calcs.
                 */
                deltzi = (deltaz + 50) / 2;

                /*
                                 * The name of this variable is lets-not-use-goto's-to-jump-out-of
                                 * -conditional-loops-whenever-we-feel-like-it;
                 */
                lnugtjooclwwfli = 0;

                do {

                    /* Check for the existence of a mixing layer. */
                    if (thick_mx == 0) {

                        /*
                                 * Calculate the partial derivatives of u and v with respect to
                                 * time when there is no mixing layer. Initialise specific humidities
                         */

 /*
                                 * If initialisation of actual wind profile and geostrophic winds do
                                 * not "match" proopery the model with crash above the height of the
                                 * mixed layer due to inertial occilation. If necessary one may take
                                 * out this feature above the height of the planetary boundary layer.
                         */
                        dudt50 = cf * (vd[1] - vgd[1])
                                - (smf
                                * (ustar * ustar * (ud[1] / awind) / 25));
                        dvdt50 = (-1) * cf * (ud[1] - ugd[1])
                                - (smf
                                * (ustar * ustar * (vd[1] / awind) / 25));
                        dv50 = dvdt50 * delta;
                        du50 = dudt50 * delta;
                        vd[1] = vd[1] + dv50;
                        ud[1] = ud[1] + du50;
                        dqdt2[1] = 0;

                        for (i = 2; i <= ntrp; i++) {
                            dudt[i] = cf * (vd[1] - vgd[i]);
                            dvdt[i] = (-1) * cf * (ud[i] - ugd[i]);
                            du_momday[i] = dudt[i] * delta;
                            dv_momday[i] = dvdt[i] * delta;
                            ud[i] = ud[i] + du_momday[i];
                            vd[i] = vd[i] + dv_momday[i];
                            dqdt2[i] = 0;
                        }

                        lnugtjooclwwfli = 1;

                        /**
                         * *******!!!!!!goto 999 *******
                         */
                    } else {

                        /* Set up contingencies to check if near top of mixing layer */
                        for (i = 2; i <= ntrp; i++) {
                            dudt[i] = cf * (vd[i] - vgd[i]);
                            dvdt[i] = (-1) * cf * (ud[i] - ugd[i]);

                            /*
                                         * if(zi[i]>hgt)
                                         * {
                                         * td[i]=td[i]-(rad-advgt)*(delta/2);
                                         * }
                             */
                            if (zi[i] <= hgt) {
                                if (zi[i] == zi[ntrp]) {
                                    km[i + 1] = 0;
                                    ud[i + 1] = ud[i];
                                    vd[i + 1] = vd[i];
                                }

                                /*
                                                 * Calc rate of change of q u and v with respect to time
                                                 * due to the diffusive properties within the layer
                                 */
                                dqdt2[i]
                                        = ((km[i + 1] * (qd[i + 1] - qd[i]) / deltaz) - (km[i] * ((qd[i] - qd[i - 1]) / deltaz)))
                                        / deltaz;
                                dudt[i]
                                        = dudt[i]
                                        + ((km[i + 1] * (ud[i + 1] - ud[i]) / deltaz) - (km[i] * (ud[i] - ud[i - 1]) / deltaz))
                                        / deltaz;
                                dvdt[i]
                                        = dvdt[i]
                                        + ((km[i + 1] * (vd[i + 1] - vd[i]) / deltaz) - (km[i] * (vd[i] - vd[i - 1]) / deltaz))
                                        / deltaz;
                            }

                            du_momday[i] = dudt[i] * (delta / 2);
                            dv_momday[i] = dvdt[i] * (delta / 2);
                            qd[i] = qd[i] + (dqdt2[i] * delta / 2);
                        }

                        /*
                                 * At lower boundary smooth out the derivatives using classicial
                                 * functional form.
                         */
                        dudt[1]
                                = cf * (vd[1] - vgd[1])
                                + (km[2] * (ud[2] - ud[1]) / (deltaz * deltzi))
                                - (smf
                                * (uif * ustar * ustar * (ud[1] / awind)
                                / deltzi));
                        dvdt[1]
                                = (-1) * cf * (ud[1] - ugd[1])
                                + (km[2] * (vd[2] - vd[1]) / (deltaz * deltzi))
                                - (smf
                                * (uif * ustar * ustar * (vd[1] / awind)
                                / deltzi));
                        dqdt2[1]
                                = (km[2] * (qd[2] - qd[1]) / (deltaz * deltzi))
                                + (evap * smf * uif) / (le * DENS * deltzi);
                        du_momday[1] = dudt[1] * (delta / 2);
                        dv_momday[1] = dvdt[1] * (delta / 2);
                        qd[1] = qd[1] + (dqdt2[1] * delta / 2);

                        /* Compute the new values for u and v components */
                        for (j = 1; j <= ntrp; j++) {
                            ud[j] = ud[j] + du_momday[j];
                            vd[j] = vd[j] + dv_momday[j];
                        }
                    }

                    if (lnugtjooclwwfli != 1) {

                        /* windspeed at level one */
                        awind = Math.sqrt(ud[1] * ud[1] + vd[1] * vd[1]);

                        /* Cycle through twice. */
 /* Another bloody goto here originally. */
                        xxx = xxx + 1;
                        over = xxx % 2;			  // was cast to (int)
                    }
                } while ((over != 0) && (lnugtjooclwwfli == 0));

                /* Cycle through twice, but only sometimes */
 /* Update the fine mesh arrays */
                /**
                 * *****************
                 * fine *****************
                 */
                _fine(u_fine, v_fine, t_fine, q_fine, ud, vd, td, qd, hgt,
                        aptemp);

                /**
                 * ********************
                 * end of momday ********************
                 */
            }

            /* Evaporative flux, surface temperature solutions */
            /**
             * ******************************************
             * flux *****************************************
             */

            /*
                         * During the day, flux calcs surface temp and surface specific humidity
                         * from temp and humidity at za and the sensible and latent heat fluxes.
                         * It also computes the updated value of the evaporative flux. During the
                         * night is calls gtemp to calc the temp in the absence of turbulence.
             */

 /*
                         * We calculate otemp from gtemp routine once we have reached
                         * radiative balance at night
             */
 /* Put in a time restraint so you do not call gtemp in the am, */
            oshum = Math.pow(10, (6.1989 - (2353 / bareradiotemp)));
            bareevapflux = le * DENS * (oshum - qd[1]) / sumw * f;

            if (qd[1] >= oshum) {
                bareevapflux = 0.001;
            }

            /**
             * *******************
             * average *******************
             */

            /* inital pass flag */
            if (average_pass == 0) {
                for (i = 1; i <= 4; i++) {
                    av_array[i] = bareevapflux;
                }

                average_pass = 1;
            } else {
                for (j = 2; j <= 4; j++) {
                    av_array[j - 1] = av_array[j];
                }

                av_array[4] = bareevapflux;
            }

            average_sum = 0;

            for (k = 1; k <= 4; k++) {
                average_sum = av_array[k] + average_sum;
            }

            evap_smooth = average_sum / 4;

            /**
             * *******************
             * end of average *******************
             */
            bareevapflux = evap_smooth;
            evap = bareevapflux;

            /* Nighttime */
            if (rnet < 0) {

                /**
                 * *******************
                 * gtemp *******************
                 */

                /*
                                 * Getemp determines the surface temperature in the adsence of
                                 * turbulence. It sets up and energy balance equation as a fourth-order
                                 * polynomial equation in surface temperature and solves it using
                                 * Newton's iteration algorithm. Note that gtemp assumes radiative
                                 * balance
                 */
                x = otemp;

                /* Set up the polynomial */
                engbal
                        = (aepsi * sigma * Math.pow((t_fine[3] - tdif_s - 1.5), 4) + swave - heat - evap)
                        * z[2] + lambda * tt[2];
                a1 = z[2] * sigma * epsi;
                a_gtemp[1] = (-1) * engbal;
                a_gtemp[2] = lambda;
                a_gtemp[3] = 0;
                a_gtemp[4] = 0;
                a_gtemp[5] = a1;

                /* Reduce the boundaries until the difference is less than 0.0001 */
                for (j = 1; j <= 20; j++) {
                    b = a_gtemp[5];
                    c = a_gtemp[5];

                    for (i = 1; i <= 3; i++) {
                        k = 5 - i;
                        b = ((0.00001 * a_gtemp[k]) + (0.00001 * x * b))
                                * 100000;
                        c = b + x * c;
                    }

                    b = a_gtemp[1] + x * b;
                    deltax = b / c;
                    x = x - deltax;
                    xdif = x - otemp;

                    /*
                                         * bad babd bad break out of the loop, but it's better than a goto
                                         * which is what was here in the fortran.
                     */
                    if (Math.abs(xdif) < 0.0001) {
                        break;
                    }
                }

                x = (otemp + x) / 2;
                otemp = x;

                /**
                 * *******************
                 * end of gtemp *******************
                 */
                bareheatflux = heat;
                vegnradiotemp = otemp;
                bareradiotemp = otemp;
                mixedradiotemp = otemp;
                evap = bareevapflux;
            } else {

                /* Daytime */
 /* bare soil */
                bareradiotemp = aptemp + (bareheatflux * sum / (DENS * CP))
                        - tdif_s;
                otemp = bareradiotemp;

                /* Daytime restraint */
                if (rnet > 0) {

                    /* Vegetation */
                    if (frveg == 1) {

                        /**
                         * **********************
                         * vegflx *********************
                         */
                        qstg = Math.pow(10, (6.1989 - 2353 / tg));
                        xleg = f * chg * DENS * le * (qstg - qaf);

                        /* xleg=(xleg+xlegn)/2; smooth not */
                        if (xleg < 0) {
                            xleg = 0;
                        }

                        tg = taf + (hg / chg) / (DENS * CP);
                        vegnevapflux = xleg + xlef;
                        taf = (cha * ta + chf * tf + tg * chg)
                                / (chf + chg + cha);
                        rprime = 1 / (1 / chf + rst);
                        qaf = (cha * qa + rprime * qstf + f * chg * qstg)
                                / (cha + rprime + f * chg);

                        /**
                         * **********************
                         * end of vegflx *********************
                         */
                        evap = vegnevapflux;
                        otemp = vegnradiotemp;

                    } else if ((frveg > 0) && (frveg < 1)) {

                        /* Mixed */
                        /**
                         * **********************
                         * vegflx *********************
                         */
                        qstg = Math.pow(10, (6.1989 - 2353 / tg));
                        xleg = f * chg * DENS * le * (qstg - qaf);

                        /* xleg=(xleg+xlegn)/2; smooth not */
                        if (xleg < 0) {
                            xleg = 0;
                        }

                        tg = taf + (hg / chg) / (DENS * CP);
                        vegnevapflux = xleg + xlef;
                        taf = (cha * ta + chf * tf + tg * chg)
                                / (chf + chg + cha);
                        rprime = 1 / (1 / chf + rst);
                        qaf = (cha * qa + rprime * qstf + f * chg * qstg)
                                / (cha + rprime + f * chg);

                        /**
                         * **********************
                         * end of vegflx *********************
                         */
                        mixedevapflux = (1 - frveg) * bareevapflux
                                + frveg * vegnevapflux;
                        mixedradiotemp
                                = Math.pow(Math.pow(bareradiotemp, 4) * (1 - frveg)
                                        + frveg * Math.pow(vegnradiotemp, 4),
                                        0.25);
                        evap = mixedevapflux;
                        otemp = mixedradiotemp;
                    }
                }
            }

            ahum = qd[1];

            /**
             * ******************************************
             * end of flux *****************************************
             */

            /* Heat flux - penman formulation */
            if (((heat >= 0) || (rnet > 0)) && (swave > 0)) {

                /**
                 * ******************************************
                 * hot *****************************************
                 */

                /*
                                 * Hot is called during the day to compute the sensible heat flux
                                 * from net all wave radiation and the evaporation
                 */
                a_sh = (lambda * (atemp - tt[2])) / z[2];
                b = lambda / (z[2] * DENS * CP);

                if ((frveg == 1) && (rnet > 0)) {

                    /* Vegetation */
                    /**
                     * **********************
                     * veghot *********************
                     */
                    hfn = DENS * CP * chf * (tf - taf);
                    hf = (hf + hfn) / 2;
                    aveg = (lambda * (taf - tt[2])) / z[2];

                    /* vgdens=ps1*100/(r*taf); */
                    hg = (rnetg - xleg - aveg) / (1 + b / chg);
                    vegnheatflux = hg + hf;

                    /**
                     * *********************
                     * co2flx *********************** Created 11/2/90f; modified
                     * 06/17/91
                     */
                    chax = ustar * ustar / (uten - uaf);
                    pex = (xlai / 2.0) + 1;
                    rair = 1.0 / chax + rzascr;
                    rroe = 1.83;
                    rafcanopy = raf * pex / xlai;

                    /* Bunch o' comments */
                    rrtot = 1.32 * rafcanopy + 1.66 * rst + rair;
                    fco2 = rroe * (co - ci) / rrtot;
                    fco2 = (fco2 * frveg) / 0.044;
                    /* In moles/m2/s */
                    ccan = co - (co - ci) * rair * frveg / rrtot;

                    /**
                     * *********************
                     * end of co2flx ***********************
                     * ********************* ozone *********************
                     */
                    chax = ustar * ustar / (uten - uaf);
                    pes = (xlai / 2.0) + 1.0;
                    rair = 1.0 / chax + rzascr;
                    rag = 1 / chg;
                    rroz = 1.9;
                    the_time = ptime;
                    rtot
                            = 1
                            / (1 / ((raf + rcut) * 1.32 / (2 * xlai))
                            + 1
                            / (raf * 1.32 * pes / xlai + rst * 1.66 * pes / xlai)
                            + 1 / (rag * 1.32));
                    caf = coz_air / (1 + rair / rtot);
                    fleaf = rroz * caf / ((raf + rcut) * 1.32 / (2 * xlai))
                            * 1.0e3;
                    fmeso = rroz * caf
                            / ((1.32 * raf * pes / xlai)
                            + (1.66 * rst * pes / xlai));
                    fmeso = fmeso * 1.0e3;
                    fg = rroz * caf / (rag * 1.32) * 1.0e3;
                    flux_plant = fleaf + fmeso;
                    fbare = rroz * coz_air / sumo3 * 1.0e3;
                    fglobal = (flux_plant + fg) * frveg + fbare * (1 - frveg);

                    /**
                     * *********************
                     * end of ozone ***********************
                     * ********************** end of veghot
                     * *********************
                     */
                    heat = vegnheatflux;
                } else if ((frveg > 0) && (frveg < 1) && (rnet > 0)) {

                    /* mixed */
                    bareheatflux = (barenetradn - bareevapflux - a_sh)
                            / (1 + b * sum);

                    /**
                     * **********************
                     * veghot 2 *********************
                     */
                    hfn = DENS * CP * chf * (tf - taf);
                    hf = (hf + hfn) / 2;
                    aveg = (lambda * (taf - tt[2])) / z[2];

                    /* vgdens=ps1*100/(r*taf); */
                    hg = (rnetg - xleg - aveg) / (1 + b / chg);
                    vegnheatflux = hg + hf;

                    /**
                     * *********************
                     * co2flx 2 *********************** Created 11/2/90f;
                     * modified 06/17/91
                     */
                    chax = ustar * ustar / (uten - uaf);
                    pex = (xlai / 2.0) + 1;
                    rair = 1.0 / chax + rzascr;
                    rroe = 1.83;
                    rafcanopy = raf * pex / xlai;

                    /* Bunch o' comments */
                    rrtot = 1.32 * rafcanopy + 1.66 * rst + rair;
                    fco2 = rroe * (co - ci) / rrtot;
                    fco2 = (fco2 * frveg) / 0.044;
                    /* In moles/m2/s */
                    ccan = co - (co - ci) * rair * frveg / rrtot;

                    /**
                     * *********************
                     * end of co2flx 2 ***********************
                     * ********************* ozone 2 *********************
                     */
                    chax = ustar * ustar / (uten - uaf);
                    pes = (xlai / 2.0) + 1.0;
                    rair = 1.0 / chax + rzascr;
                    rag = 1 / chg;
                    rroz = 1.9;
                    the_time = ptime;
                    rtot
                            = 1
                            / (1 / ((raf + rcut) * 1.32 / (2 * xlai))
                            + 1
                            / (raf * 1.32 * pes / xlai + rst * 1.66 * pes / xlai)
                            + 1 / (rag * 1.32));
                    caf = coz_air / (1 + rair / rtot);
                    fleaf = rroz * caf / ((raf + rcut) * 1.32 / (2 * xlai))
                            * 1.0e3;
                    fmeso = rroz * caf
                            / ((1.32 * raf * pes / xlai)
                            + (1.66 * rst * pes / xlai));
                    fmeso = fmeso * 1.0e3;
                    fg = rroz * caf / (rag * 1.32) * 1.0e3;
                    flux_plant = fleaf + fmeso;
                    fbare = rroz * coz_air / sumo3 * 1.0e3;
                    fglobal = (flux_plant + fg) * frveg + fbare * (1 - frveg);

                    /**
                     * *********************
                     * end of ozone 2 ***********************
                     * ********************** end of veghot 2
                     * *********************
                     */
                    mixedheatflux = (1 - frveg) * bareheatflux
                            + frveg * vegnheatflux;
                    heat = mixedheatflux;
                } else {

                    /* bare soil */
                    bareheatflux = (barenetradn - bareevapflux - a_sh)
                            / (1 + b * sum);
                    heat = bareheatflux;
                }

                /**
                 * ******************************************
                 * end of hot *****************************************
                 */
            }

            if (rnet > 0) {
                stabcriteria = 1;
            } else {
                stabcriteria = 2;
            }

            /* End of the atmospheric cycle */
            /**
             * ******************************************
             * below *****************************************
             */

            /*
                         * Below is called every time step (N+1) to update the sub-surface
                         * temps. The lowest level has a constant time step btemp, whereas the
                         * surface temp has been found previously in flux. Use of the
                         * lea-frog method computes sub-surface temps at the next time
                         * step from those at the current and the previous time step. below
                         * also calls water to update the sub-surface soil moisture status.
             */
 /* tt[2] is the tmperature at the first level below the soil */
            nlvl1 = nlvls + 1;

            /*
                         * Use the fraction of vegetation to set the boundary conditions
                         * for the first level in the ground (surface).
             */
            if (frveg == 0) {
                tt[1] = otemp;
            } else if ((frveg > 0) && (frveg < 1)) {
                tt[1] = Math.pow(Math.pow(tg, 4) * frveg
                        + Math.pow(bareradiotemp, 4) * (1 - frveg),
                        0.25);
            } else if (frveg == 1) {
                tt[1] = tg;
            } else {

                /* Some kind of error message here. */
                // System.out.println("Error in below - bad value for frveg." + frveg + "\n");
                throw new Exception("Error in below - bad value for frveg. %f\n"
                        + frveg);
            }

            if ((heat < 0) || (rnet < 0)) {
                tt[1] = otemp;
            }

            /* k refers to level, k=1 being the surface */
            if (time == 0) {
                for (i = 1; i <= nlvl1; i++) {
                    te[i] = tt[i];
                }
            }

            /*
                         * tt is at time N, te is at time n-1, ttt is at time n+1.
                         * Here we represent the diffusion term. See non-existent manual
                         * explanation
             */
            for (k = 2; k <= nlvls; k++) {
                term1 = kappa / (xfun[k] * xfun[k] * del * del);
                term2 = (tt[k + 1] - 2 * te[k] + tt[k - 1]) / (dzeta * dzeta);
                term3 = (te[k + 1] - tt[k - 1]) / 2 * dzeta;
                dtdt[k] = term1 * (term2 - term3);

                if (time == 0) {

                    /* Initial computation of tt[k] and redefinition for the future */
                    ttt[k] = te[k] + delta * dtdt[k];
                    tt[k] = ttt[k];
                } else {

                    /*
                                         * Computation of ttt[k] via leapfrog rule, along with redefinition
                                         * for the future
                     */
                    ttt[k] = te[k] + 2 * delta * dtdt[k];
                    te[k] = tt[k];
                    tt[k] = ttt[k];
                }
            }

            /* I have no clue why this is here */
            if (ptime == 16) {
                dummy = 1;
            }

            /* Skip the substrate water component is wmax>1 */
            if (wmax > 1) {
                w2g = 99.99;
                wgg = 99.99;
            } else {

                /**
                 * *********************
                 * water *********************
                 */

                /*
                                 * Water is baed on the technique of Deardorff (1978.) It uses the
                                 * evaporative flux value obtained in flux and updates two internal
                                 * variables wgg and w2g, which represent the soil moisture content
                                 * of the soil close to the surface and in the first 50cm of soil,
                                 * respectively. The empirical constants can be found in the article.
                 */

 /*
                                 * 1988 modifications: substrate layer assigned moisture availability
                                 * which is called fsub (w2g/wmax). Intermediate layer water equation
                                 * for variable called (win). Top layer is assumed to pertain to top
                                 * 2cm (instead of top 10cm as for intermediate layer), 20% of root
                                 * evaporation is drawn from this layer. Top layer draws on all surface
                                 * evaporation. Lowest (reservoir) layer on all evaporation. Note that
                                 * constsnt const2 changed from earlier version.
                 */

 /*
                                 * Note if you wish to suppress variation in substrate water content
                                 * with time, let wmax equal to a very large value, e.g. 10
                 */
                if (time == 0) {
                    win = (wgg + w2g) / 2;
                }

                per = OMG * 3600;
                c11 = CONST1 / (RHOW * DIP);
                c22 = CONST2 / per;
                c33 = 1 / (D2P * RHOW);
                c44 = CONST1 / (RHOW * DINT);
                evax = evap / le;

                if ((frveg > 0) && (rnet > 0)) {
                    evas = (xleg * frveg + (1 - frveg) * bareevapflux) / le;
                    evai = (xlef * frveg) / le;
                } else {
                    evas = evax;
                    evai = 0;
                }

                ww1 = (c11 * evas + c22 * (wgg - win)) * delta;
                ww2 = (c44 * evai - c22 * (wgg + w2g - 2 * win)) * delta;
                ww3 = evai * c33 * delta;
                wgg = wgg - ww1;
                win = win - ww2;
                w2g = w2g - ww3;

                if (wgg <= 0) {
                    wgg = 0.001;
                }

                if (win <= 0) {
                    win = 0.001;
                }

                if (w2g <= 0) {
                    w2g = 0.001;
                }

                /*
                                 * compute the updated version of moisture availability and substrate
                                 * moisture availability
                 */
                f = (wgg / wmax); //ok let's NOT check this out...
                fsub = (w2g / wmax);

                /**
                 * *********************
                 * end of water *********************
                 */
            }

            /**
             * ******************************************
             * end of below *****************************************
             */

            /* output is written every outtt seconds */
 /* if(tmod==0) */
            // For debug use next line
            // if (out_count < output_data[0].length)
            if (tmod == 0) {
                /**
                 * update progress monitor in main program
                 */
                /*update by Yannis Konstas(start)*/
                if (!batch) {
                    //psubamsRef.updateModelProgress();
                }
                /*update by Yannis Konstas(end)*/
                /**
                 * ******************************************
                 * output *****************************************
                 */

                /* Finally output some variables */
 /* Initialise some variables in the usual brain-dead way */
                if (output_pass == 0) {
                    i_header_output = 1;
                    undefined = Double.NaN;
                    output_pass = 1;
                }
                undefined = Double.NaN;
                g_flux = rnet - heat - evap;
                bowen = heat / evap;

                if (bowen < 0.0) {
                    bowen = undefined;
                }

                if (frveg != 0) {
                    air_leaf_t = taf - 273.23;

                    /* ground_t=tg-273.23 */
                } else {
                    air_leaf_t = undefined;

                    /* ground_t=undefined */
                    vfl = undefined;
                }

                pes = (xlai / 2.0) + 1;
                stom_r = rs * pes / xlai;
                co2_flux = fco2 * 1.0e6;
                ccan_concn = ccan * 1.0e6;
                water_use_eff = (co2_flux * 4.4e-8) / (xlef / le);

                // To change the variables that are undefined at the end of daytime profile.
                if (rnet < 0) {

                    stom_r = undefined;
                    ccan_concn = undefined;
                    water_use_eff = undefined;
                    caf = undefined;
                    flux_plant = undefined;
                    vfl = undefined;
                    uaf = undefined;
                    co2_flux = undefined;
                    psie = undefined;
                    psim = undefined;
                    fglobal = undefined;
                    //swave = undefined;
                }

                /*update by Yannis Konstas(start)*/
                //calculate EF1, EF2 and VMC (Volumetric Moisture Content)
                //evap = 0.4d; heat = 0.7d;
                ef1 = evap / (evap + heat);
                ef2 = heat / (heat + evap);

                //calculate Ed (LEdaily)
                double edTime = realtm * (1.0 / (86400d));
                ed1 = (2 * edN / (Math.PI * Math.sin((Math.PI * edTime) / edN))) * evap;
                ed1 *= 0.0864; // W*m-2 -> MJm-2d-1
                /*update by Yannis Konstas(end)*/
 /*
                                 * TCS moved names of output data to RunModel
                                 * if(i_header_output==1)
                                 * {
                                 * i_columns=OUTPUT_COUNT;
                 */
 /* Print out the headers */

 /*
                                 * fprintf(o_model,"%d %d\n",i_columns,no_rows);
                                 * fprintf(o_model,"time\n");
                                 * fprintf(o_model,"shortwave_flux/Wm-2\n");
                                 * fprintf(o_model,"net_radiation/Wm-2\n");
                                 * fprintf(o_model,"sensible_heat_flux/Wm-2\n");
                                 * fprintf(o_model,"latent_heat_flux/Wm-2\n");
                                 * fprintf(o_model,"ground_flux/Wm-2\n");
                                 * fprintf(o_model,"air_temperature_50m/C\n");
                                 * fprintf(o_model,"air_temperature_10m/C\n");
                                 * fprintf(o_model,"air_temperature_foliage/C\n");
                                 * fprintf(o_model,"radiometric_temperature/C\n");
                                 * fprintf(o_model,"wind_50m/kts\n");
                                 * fprintf(o_model,"wind_10m/kts\n");
                                 * fprintf(o_model,"wind_in_foliage/kts\n");
                                 * fprintf(o_model,"specific_humidity_50m/gKg-1\n");
                                 * fprintf(o_model,"specific_humidity_10m/gKg-1\n");
                                 * fprintf(o_model,"specific_humidity_foliage/gKg-1\n");
                                 * fprintf(o_model,"bowen_ratio\n");
                                 * fprintf(o_model,"surface_moisture_availability\n");
                                 * fprintf(o_model,"root_zone_moisture_availability\n");
                                 * fprintf(o_model,"stomatal_resistance/sm-1\n");
                                 * fprintf(o_model,"vapor_pressure_deficit/mbar\n");
                                 * fprintf(o_model,"leaf_water_potential/bars\n");
                                 * fprintf(o_model,"epidermal_water_potential/bars\n");
                                 * fprintf(o_model,"ground_water_potential/bars\n");
                                 * fprintf(o_model,"co2_flux/micromolesm-2s-1\n");
                                 * fprintf(o_model,"co2_concentration_canopy/ppmv\n");
                                 * fprintf(o_model,"water_use_efficiency\n");
                                 * fprintf(o_model,"o3_concentration_canopy/ppmv\n");
                                 * fprintf(o_model,"global_o3_flux/ugm-2s-1\n");
                                 * fprintf(o_model,"o3_flux_plant/ugm-2s-1\n");
                
                 */
 /* Copy the output names to the transfer buffers */
 /*
                                 * output_names[0] = OUT_TIME;
                                 * output_names[1] = OUT_SHORTWAVE_FLUX;
                                 * output_names[2] = OUT_NET_RADIATION;
                                 * output_names[3] = OUT_SENSIBLE_HEAT_FLUX;
                                 * output_names[4] = OUT_LATENT_HEAT_FLUX;
                                 * output_names[5] = OUT_GROUND_FLUX;
                                 * output_names[6] = OUT_AIR_TEMP_50M;
                                 * output_names[7] = OUT_AIR_TEMP_10M;
                                 * output_names[8] = OUT_AIR_TEMP_FOLIAGE;
                                 * output_names[9] = OUT_RADIOMETRIC_TEMP;
                                 * output_names[10] = OUT_WIND_50M;
                                 * output_names[11] = OUT_WIND_10M;
                                 * output_names[12] = OUT_WIND_IN_FOLIAGE;
                                 * output_names[13] = OUT_SPEC_HUMID_50M;
                                 * output_names[14] = OUT_SPEC_HUMID_10M;
                                 * output_names[15] = OUT_SPEC_HUMID_FOLIAGE;
                                 * output_names[16] = OUT_BOWEN_RATIO;
                                 * output_names[17] = OUT_SURF_MOIST_AVAIL;
                                 * output_names[18] = OUT_ROOT_ZONE_MOIST_AVAIL;
                                 * output_names[19] = OUT_STOMATAL_RESIST;
                                 * output_names[20] = OUT_VAPOR_PRESS_DEFICIT;
                                 * output_names[21] = OUT_LEAF_WATER_POTENTIAL;
                                 * output_names[22] = OUT_EPIDERMAL_WATER_POT;
                                 * output_names[23] = OUT_GROUND_WATER_POT;
                                 * output_names[24] = OUT_CO2_FLUX;
                                 * output_names[25] = OUT_CO2_CONC_CANOPY;
                                 * output_names[26] = OUT_WATER_USE_EFF;
                                 * output_names[27] = OUT_O3_CONC_CANOPY;
                                 * output_names[28] = OUT_GLOBAL_O3_FLUX;
                                 * output_names[29] = OUT_O3_FLUX_PLANT;
                
                                 * output_names[30] = OUT_WET_BLB_POT_TEMP;
                                 * output_names[31] = OUT_LFT_INDEX;
                                 * output_names[32] = OUT_K_INDEX;
                                 * output_names[33] = OUT_TOT_TOTLS;
                                 * output_names[34] = OUT_SWT_INDEX;
                                 * 
                                 * 
                                 * i_header_output=2;
                                 * }
                 */
 /* Print out the output data */
 /*
                                 * fprintf(o_model,"%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.1f %.2f %.2f %.2f %.3f %.2f %.2f %.3f %.3f %.3f %.3f\n",
                                 * ptime,swave,rnet,heat,evap,g_flux,atemp-273.23,ta-273.23,air_leaf_t,otemp-273.23,
                                 * awind*1.98,uten*1.98,uaf*1.98,q_fine[1]*1000,qa*1000,qaf*1000,bowen,
                                 * f,fsub,stom_r,vfl,psim,psie,psig,co2_flux,ccan_concn,water_use_eff,
                                 * caf,fglobal,flux_plant);
                 */
 /* Copy the output data into the transfer buffer */
                output_data[0][out_count] = ptime;
                output_data[1][out_count] = swave;
                output_data[2][out_count] = rnet;
                output_data[3][out_count] = heat;
                output_data[4][out_count] = evap;
                output_data[5][out_count] = g_flux;
                output_data[6][out_count] = atemp - 273.23;
                output_data[7][out_count] = ta - 273.23;
                output_data[8][out_count] = air_leaf_t;
                output_data[9][out_count] = otemp - 273.23;
                output_data[10][out_count] = awind * 1.98;
                output_data[11][out_count] = uten * 1.98;
                output_data[12][out_count] = uaf * 1.98;
                output_data[13][out_count] = q_fine[1] * 1000;
                output_data[14][out_count] = qa * 1000;
                output_data[15][out_count] = qaf * 1000;
                output_data[16][out_count] = bowen * 1000;
                output_data[17][out_count] = f * 1000;
                output_data[18][out_count] = fsub * 1000;
                output_data[19][out_count] = stom_r;
                output_data[20][out_count] = vfl;
                output_data[21][out_count] = psim;
                output_data[22][out_count] = psie;
                output_data[23][out_count] = psig * 1000;
                output_data[24][out_count] = co2_flux;
                output_data[25][out_count] = ccan_concn;
                output_data[26][out_count] = water_use_eff * 1000;
                output_data[27][out_count] = caf * 1000;
                output_data[28][out_count] = fglobal;
                output_data[29][out_count] = flux_plant;

                /* And do the same for the sounding data */
 /* Do it for all levels; this is every 50 meters */
                for (i = 0; i < SOUNDING_LEVELS; i++) {

                    /* Height */
                    sounding_data[0][out_count][i] = 50 * i;

                    /* Pressure */
                    sounding_data[1][out_count][i] = ps[1] * Math.exp(-9.8 * 50 * i / (287 * (ts[1] + C_TO_K)));

                    /* Temperature */
                    sounding_data[2][out_count][i] = t_fine[i + 1] - 273.15;

                    /* Humidity */
                    sounding_data[3][out_count][i] = q_fine[i + 1] * 1000;

                    /* u component of wind */
                    sounding_data[4][out_count][i] = u_fine[i + 1] * 1.98;

                    /* v component of wind */
                    sounding_data[5][out_count][i] = v_fine[i + 1] * 1.98;

                    /* Potential temperature */
                    sounding_data[6][out_count][i] = t_fine[i + 1] * Math.pow((1000 / sounding_data[1][out_count][i]), 0.286) - 273.15;

                    /* Magnitude of wind */
                    sounding_data[7][out_count][i] = Math.sqrt(Math.pow(u_fine[i + 1], 2) + Math.pow(v_fine[i + 1], 2));

                    /* Direction of wind */
                    sounding_data[8][out_count][i] = (Math.acos(u_fine[i + 1] / sounding_data[7][out_count][i])) * RADIANS_TO_DEGREES;

                }


                /*	Calculations for New Indices */
 /* Calculate some new indices. This was supposed to be trivial, only
                                it isn't. The new indices are mostly calculated in terms of temperatures
                                (and/or dew point potential temperatures) at a given pressure
                                level. This would be okay if we could use the x_fine[] value, but
                                these don't go up very high, so we have to resort to using the
                                regular non-interpolated sounding data as well, combining the two,
                                and tabulating that against the pressure. */
 /* Get the wet bulb potential temperature */
 /* Calculate the averages of the dew point and temperature at the
                                two lowest levels first, because we're going to need them in a bit. */
                dp_average = (_dpt_from_qp(q_fine[1], sounding_data[1][out_count][0]) + _dpt_from_qp(q_fine[2], sounding_data[1][out_count][1])) / 2.0;
                t_average = (sounding_data[2][out_count][0] + sounding_data[2][out_count][1]) / 2 + C_TO_K;
                wbpt = _wetblb((sounding_data[1][out_count][0] + sounding_data[1][out_count][1]) / 2, t_average, dp_average);

                //output_data[30][out_count]= wbpt;
                /* Now the tricky bastards. The arrays we're going to use are:
                                ud[], vd[], td[], qd[] amd u_fine[], v_fine[], t_fine[], q_fine[]
                                u,v,t,q and q refer to:
                                u -> u component of wind in ms-1
                                v -> v component of wind in ms-1
                                t -> temperature in kelvin.
                                q -> humidity in g/g.
                                ud, vd, td, and qd start at 50m, and are every 250m from there on up.
                                u_fine, v_fine, t_fine, and q_fine start at 50m and at every 50m from there on up.
                                There are always 51 values for t_fine &c, unit offset indexing.
                                The number of ud[] values depends on the sounding - it's given by ntrp. */
 /* get the showalter index, except it isn't really, it's the "lifted index".
                                calculate some stuff that we can re-use first. */
                tat500 = _val_at_pressure(500, t_fine, td, SOUNDING_LEVELS, ntrp, ts[1], ps[1]) - C_TO_K;

                //output_data[31][out_count]= _ps_show(t_average - C_TO_K, dp_average - C_TO_K, tat500);           
                /* Functions to calculate the showalter index. */
 /* This function calculates the showalter index.
                                show=t500 - tparcel
                                tparcel=Temperature in celcius at 500mb of a parcel 
                                which lies on the moist adiabat determined by the
                                sounding at 850mb.
                 */
                p850 = 850;

                /* Find equivalent potential temperature at the LCL
                                   using 850mb temperature and dewpoint */
                thtlcl = _pr_thte(p850, t_average - C_TO_K, dp_average - C_TO_K);

                p500 = 500.0;

                /* Find the parcel temperature along pseudoadiabat at 500mb */
                guess = 0;

                /* Ah shit where did the p500 come from? 
                                   Figured it out. p500 is the pressure at 500mb... which is sort
                                   of obliged to be 500mb. God only knows why they did this. */
                pr_tmst = _pr_tmst(thtlcl, p500, guess);

                /* Subtract the parcel temp. from the temp at 500mb. 
                 */
                //t500k = (tat500 + C_TO_K );
                ps_show = (tat500 + C_TO_K) - pr_tmst;

                /* Next one, the "k index." */
                tat850 = _val_at_pressure(850, t_fine, td, SOUNDING_LEVELS, ntrp, ts[1], ps[1]) - C_TO_K;
                tat700 = _val_at_pressure(700, t_fine, td, SOUNDING_LEVELS, ntrp, ts[1], ps[1]) - C_TO_K;
                tdat850 = _dpt_from_qp(_val_at_pressure(850, q_fine, qd, SOUNDING_LEVELS, ntrp, ts[1], ps[1]), 850) - C_TO_K;
                tdat700 = _dpt_from_qp(_val_at_pressure(700, q_fine, qd, SOUNDING_LEVELS, ntrp, ts[1], ps[1]), 700) - C_TO_K;

                //output_data[32][out_count] = tat850 + tdat850 + tdat700 - tat700 - tat500;
                /* Now the "total totals" */
                //output_data[33][out_count] = tat850 + tdat850 - 2 * tat500;
                /* Finally, the "sweat index"
                                Get some of the stuff that we need - the winds at 850 and 500 mb, in kts */
                uat850 = _val_at_pressure(850, u_fine, ud, SOUNDING_LEVELS, ntrp, ts[1], ps[1]);
                vat850 = _val_at_pressure(850, v_fine, vd, SOUNDING_LEVELS, ntrp, ts[1], ps[1]);
                uat500 = _val_at_pressure(500, u_fine, ud, SOUNDING_LEVELS, ntrp, ts[1], ps[1]);
                vat500 = _val_at_pressure(500, v_fine, vd, SOUNDING_LEVELS, ntrp, ts[1], ps[1]);
                wsat850 = Math.sqrt(uat850 * uat850 + vat850 + vat850);
                wsat500 = Math.sqrt(uat500 * uat500 + vat500 * vat500);

                /* Now figure out the wind direction, in degrees */
 /* Note that atan2 (in C) is atan2(_y_,_x_) which will return a value 
                                between pi and -pi, measuring anticlockwise from the x-axis, like 
                                mathematicians do. If one switches the co-ordinates (like I'm going to,
                                atan2(u,v) then the answer is conveniently clockwise from north, like
                                meteorologists do. */
                wdir850 = 180 * Math.atan2(uat850, vat850) / PI;
                wdir500 = 180 * Math.atan2(uat500, vat500) / PI;

                /* Make them between 0 and 360 */
                if (wdir850 < 0) {
                    wdir850 = 360 + wdir850;
                }
                if (wdir500 < 0) {
                    wdir500 = 360 + wdir500;
                }

                /* Calculate the completely retarded correction factor if a set of
                                conditions are met, otherwise the correction is just zero. */
                swcorr = 0;
                if ((wdir850 >= 130) && (wdir850 <= 250) && (wdir500 >= 210) && (wdir500 <= 310) && ((wdir500 - wdir850) >= 0) && ((wdir850 + wdir500) >= 15)) {
                    swcorr = Math.sin(PI * (wdir500 - wdir850) / 180.0);
                }

                //output_data[34][out_count] = 12 * tdat850 + 20 * (output_data[33][out_count] - 49.0) + 2 * wsat850 + wsat500 + 125 * swcorr;
                // New Indices.
                output_data[30][out_count] = wbpt;
                output_data[31][out_count] = ps_show;
                output_data[32][out_count] = tat850 + tdat850 + tdat700 - tat700 - tat500;
                output_data[33][out_count] = tat850 + tdat850 - 2 * tat500;
                output_data[34][out_count] = 12 * tdat850 + 20 * ((tat850 + tdat850 - 2 * tat500) - 49.0)
                        + 2 * wsat850 + wsat500 + 125 * swcorr;	// used o/p data[33] earlier

                /*update by Yannis Konstas(start)*/
                output_data[35][out_count] = ef1;
                output_data[36][out_count] = ef2;

                //if (evap + heat > 0)
                if (rnet - g_flux > 0 && ef1 < 3 && ef2 < 3) {
                    data_graphs[0][out_count] = ef1;
                    data_graphs[1][out_count] = ef2;
                } else {
                    data_graphs[0][out_count] = Double.NaN;
                    data_graphs[1][out_count] = Double.NaN;
                }
                output_data[37][out_count] = lwdn;

                /*update by Yannis Konstas(start)*/
 /* uplong */
                if (frveg == 1.0) {
                    if (rnet > 0) {
                        lwup = epsi * sigma * Math.pow(vegnradiotemp, 4);
                    } else {
                        lwup = undefined;
                    }
                }
                /*update by Yannis Konstas(end)*/

                output_data[38][out_count] = lwup;
                output_data[39][out_count] = ed1;

                if (ed1 > 0 && ed1 < 40) //rnet
                {
                    data_graphs[2][out_count] = ed1;
                } else {
                    data_graphs[2][out_count] = Double.NaN;
                }

                /*update by Yannis Konstas(end)*/
                out_count += 1;

                /**
                 * ******************************************
                 * end of output *****************************************
                 */
            }

            /* Increment time */
            time += delta;
        } while (realtm < timend);
        /* End of the "do" statement */


        // To Remove the first calculation and remove the initial distortions in graphs. 06/02/03 : sst
        for (int params = 0; params < OUTPUT_COUNT; params++) {
            output_data[params][0] = output_data[params][1];	// For output data array
        }

        for (int sndgs = 0; sndgs < SOUNDING_COUNT; sndgs++) {
            for (int sndlvl = 0; sndlvl < SOUNDING_LEVELS; sndlvl++) {		// For sounding array
                sounding_data[sndgs][0][sndlvl] = sounding_data[sndgs][1][sndlvl];
            }
        }


        /* Close the output file */
        // tcs fclose(o_model);

        /* done */
        // return(0);
    }

    // Functions begin here

    /*
	 * Returns the critical stomatal resistance for the critical
	 * ground water potential.
     */
    /**
     * Method declaration
     *
     *
     * @param b1
     * @param psice
     * @param rmin
     * @param fs
     * @param ft
     *
     * @return
     *
     * @see
     */
    public double _stomc(double b1, double psice, double rmin, double fs,
            double ft) {
        double fpsice;

        fpsice = 1 + b1 * psice;

        return (rmin * fs * fpsice * ft);
    }

    /* The exponential function for solar radiation - Albert Olioso. */
    /**
     * Method declaration
     *
     *
     * @param sc
     * @param sol
     *
     * @return
     *
     * @see
     */
    public double _stomfs(double sc, double sol) {
        return ((1 / (1 - Math.exp(-1 / sc * sol))));
    }

    /*
	 * Calculate solar transmissions using the look up tables produced
	 * in gettbl.
	 * 
	 * Originally this was a single function, now split into two (one of
	 * which will be called twice - with different parameters each time)
	 * so no more duplication of functionality, and no more global variables.
     */

 /*
	 * This one can be used for either ftabs or ftscat, since the
	 * same things are done to both of them.
     */
    /**
     * Method declaration
     *
     *
     * @param path
     * @param tbl
     * @param psi
     *
     * @return
     *
     * @see
     */
    public double _transmft(double path, double[] tbl, double psi) {
        double f;

        if (path >= 10.0) {
            f = tbl[46];
        } else {
            f = _transmint(tbl, path);
            f = (psi / 1013.25) * (f - 1) + 1;
        }

        return (f);
    }

    /**
     * Method declaration
     *
     *
     * @param path
     * @param tbl
     *
     * @return
     *
     * @see
     */
    public double _transmfb(double path, double[] tbl) {
        double f;

        if (path >= 10.0) {
            f = tbl[46];
        } else {
            f = _transmint(tbl, path);
        }

        return (f);
    }

    /**
     * Method declaration
     *
     *
     * @param table
     * @param path
     *
     * @return
     *
     * @see
     */
    public double _transmint(double[] table, double path) {
        int ipath;
        double fract;

        fract = 5 * (path - 1) + 1;
        ipath = (int) fract;
        fract = fract - ipath;

        return (_interpol(table[ipath], table[ipath + 1], fract));
    }

    /**
     * Method declaration
     *
     *
     * @param a
     * @param b
     * @param f
     *
     * @return
     *
     * @see
     */
    public double _interpol(double a, double b, double f) {
        return ((f * b) + ((1 - f) * a));
    }

    /**
     * Method declaration
     *
     *
     * @param rks
     * @param thv
     * @param thmax
     * @param cosbyb
     *
     * @return
     *
     * @see
     */
    public double _cond(double rks, double thv, double thmax, double cosbyb) {
        double rkw;

        rkw = ((6.9e-6) * rks
                * Math.pow((thv / (thmax * 0.75)), (2 * cosbyb + 2)));

        return (rkw);
    }

    /* Function definitions used in vel. */
    /**
     * Method declaration
     *
     *
     * @param wind
     * @param height
     * @param roughness
     * @param stability
     *
     * @return
     *
     * @see
     */
    public double _you_star(double wind, double height, double roughness,
            double stability) {
        return ((KARMAN * wind
                / ((Math.log(height / roughness) + stability))));
    }

    /**
     * Method declaration
     *
     *
     * @param friction
     * @param height
     * @param roughness
     * @param stability
     *
     * @return
     *
     * @see
     */
    public double _r_ohms(double friction, double height, double roughness,
            double stability) {
        return ((0.74 * (Math.log(height / roughness) + stability)
                / (KARMAN * friction)));
    }

    /**
     * Method declaration
     *
     *
     * @param star
     * @param height
     * @param roughness
     * @param stability
     *
     * @return
     *
     * @see
     */
    public double _wind(double star, double height, double roughness,
            double stability) {
        return (((star / KARMAN)
                * (Math.log(height / roughness) + stability)));
    }

    /**
     * Method declaration
     *
     *
     * @param height
     * @param mol
     *
     * @return
     *
     * @see
     */
    public double _stab(double height, double mol) {
        return ((Math.pow((1 - 15 * height / mol), 0.25)));
    }

    /**
     * Method declaration
     *
     *
     * @param height
     * @param mol
     *
     * @return
     *
     * @see
     */
    public double _stabh(double height, double mol) {
        return ((Math.pow((1 - 9 * height / mol), 0.5)));
    }

    /**
     * Method declaration
     *
     *
     * @param par1
     * @param par2
     *
     * @return
     *
     * @see
     */
    public double _fstabh(double par1, double par2) {
        return ((2 * Math.log((par1 + 1) / (par2 + 1))));
    }

    /**
     * Method declaration
     *
     *
     * @param par1
     * @param par2
     *
     * @return
     *
     * @see
     */
    public double _fstabm(double par1, double par2) {
        return ((Math.log(((par1 * par1 + 1) * Math.pow((par1 + 1), 2)) / ((par2 + 1) * Math.pow((par2 + 1), 2)))
                + 2 * (Math.atan(par2) - Math.atan(par1))));
    }

    /**
     * Method declaration
     *
     *
     * @param star
     * @param roughness
     * @param par3
     *
     * @return
     *
     * @see
     */
    public double _restrn(double star, double roughness, double par3) {
        return (((Math.log(KARMAN * star * roughness + par3) - Math.log(par3))
                / (KARMAN * star)));
    }

    /*
	 * Function to computer eddy diffusivities as a function of height,
	 * friction velocity and Monin Obukhov length for all the relevant levels
	 * in the mixing layer at each time step during the day; using the method
	 * of O'Brien. Used in Momday.
     */
    /**
     * Method declaration
     *
     *
     * @param km
     * @param zi
     * @param zk
     * @param hgt
     * @param za
     * @param ustar
     * @param mol
     * @param ntrp
     * @param deltaz
     *
     * @return
     *
     * @see
     */
    public double _daykm(double[] km, double[] zi, double[] zk, double hgt,
            double za, double ustar, double mol, int ntrp,
            double deltaz) {
        double ktop, kma, kmapri, thick;

        /*
		 * Array declared of size 51 so we can use fortran-style
		 * unit offset indexing.
         */
        double[] kmw = new double[51];
        int m, i;

        thick = hgt - za;

        /* only do this if there is a mixing layer. */
        if (thick != 0) {

            /*
	 * Define or calc eddy diff's (K's) at the top of the
	 * surface layer using the standard flux profile law.
	 * Calculate the derivative of kma and define k at the top
	 * of the mixing layer to be zero.
             */
            kma = (KARMAN * ustar * za
                    * Math.pow((1.0 - (15 * za / mol)), 0.25));
            km[1] = kma;
            kmw[1] = kma;
            kmapri = (kma
                    * (1 / za - (15.0 / mol) / (1.0 - 15.0 * za / mol)));
            ktop = 0;

            /*
	 * Calc. heights between the 250 metre intervals (zi
	 * system)	as passed from spline, and call zk system
	 * (interval 50,175,425 etc.)
             */
            m = ntrp + 1;

            for (i = 2; i <= m; i++) {
                zk[i] = (zi[i] - 0.5 * deltaz);
            }

            zk[1] = zi[1];

            /*
	 * Calc. the eddy diffusivities for momentum and water
	 * up to the top of the mixing layer using the O'Brien
	 * function. If ntrp goes abovehgt, fill with zeros.
             */
            for (i = 2; i <= ntrp; i++) {
                if (zk[i] < hgt) {
                    km[i] = (ktop
                            + (Math.pow((zk[i] - hgt), 2) / (thick * thick))
                            * (kma - ktop
                            + (zk[i] - za)
                            * (kmapri + 2 * (kma - ktop) / thick)));
                    kmw[i] = (ktop
                            + (Math.pow((zi[i] - hgt), 2) / (thick * thick))
                            * (kma - ktop
                            + (zi[i] - za)
                            * (kmapri + 2 * (kma - ktop) / thick)));
                } else {
                    km[i] = 0;
                    kmw[i] = 0;
                }
            }

            /* Smooth k's as a weighted average */
            for (i = 2; i <= ntrp; i++) {
                km[i] = ((kmw[i] + km[i] + kmw[i - 1]) / 3.0);
            }
        }
        /* End of thick!=0 conditional */

        return (thick);
    }

    /**
     * Method declaration
     *
     * Converts a time given in the completely stupid and annoying hour.minute
     * format into a proper decimal, in a completely stupid and annoying
     * non-obvious manner. I really really hate implicit declarations.
     *
     *
     * @param timin
     *
     * @return
     *
     * @see
     */
    public double _dectim(double timin) {
        int intim;
        double rintim;

        /*
		 * Zero ut the two least significant integral digits and any fractional
		 * component of timin
         */
        intim = (int) (timin / 100);
        rintim = (intim * 100.0);

        /* Returns the minutes part of time*60+ hour part of time*3600 */
        return ((((timin - rintim) / 60.0 + intim) * 3600.0));
    }

    /**
     * Method declaration
     *
     *
     * @param a
     * @param b
     *
     * @return
     *
     * @see
     */
    public double _min(double a, double b) {
        return ((a < b) ? a : b);
    }

    public double _max(double x, double y) {
        return ((x > y) ? x : y);
    }

    /*	Methods introduced for the new indices.	*/

 /* 	Get the value of some sounding variable at a particular pressure level from that variable at given heights, from
	*	the regular and the fine versions of the sounding. We're going to do this just by linearly interpolating between
	*	the two closest points, which possibly isn't the best way but it should be close enough, especially since this is
	*	only used for some vague and probably meaningless indexes that no-one cares about. (Which begs the question as to
	*	why they're being included at all. 
     */
    public double _val_at_pressure(double pres, double y_fine[], double y_coarse[], int yf_count, int yc_count, double base_t, double base_p) {

        double retv = 0.0;
        double p1, p2;
        int i;
        int flag = 0;

        /* NOTE: We need to flag to check if retv has been calculated, because
		it's value could be just about anything. */

 /* Start out with the fine array, see if we can get the pressure low enough;
		If the first pressure level is already lower than pres, then just return
		the first value from the fine array.*/
        if (base_p < pres) {
            retv = y_fine[1];
            flag = 1;
        }
        /* Now try the rest of the array. See if we can find pressure levels that
		bound the one we're looking for. The fine arrays will always have values
		for indexes from 1-51, but we'll use SOUNDING_LEVELS just in case */

        for (i = 0; i < yf_count - 1; i++) {
            p1 = base_p * Math.exp(-9.8 * 50 * i / (287 * (base_t + C_TO_K)));
            p2 = base_p * Math.exp(-9.8 * 50 * (i + 1) / (287 * (base_t + C_TO_K)));
            if ((p1 >= pres) && (p2 < pres)) {

                /* Got one. Interpolate. */
                retv = y_fine[i] + (y_fine[i + 1] - y_fine[i]) * (pres - p1) / (p2 - p1);
                flag = 1;
            }
        }

        /* See if we have a value already; if we don't, then do the same thing
		but with the coarse array */
        if (flag == 0) {
            for (i = 0; i < yc_count - 1; i++) {
                p1 = base_p * Math.exp(-9.8 * (50 + 250 * i) / (287 * (base_t + C_TO_K)));
                p2 = base_p * Math.exp(-9.8 * (50 + 250 * (i + 1)) / (287 * (base_t + C_TO_K)));
                if ((p1 >= pres) && (p2 < pres)) {

                    /* Got one. Interpolate */
                    retv = y_coarse[i] + (y_coarse[i + 1] - y_coarse[i]) * (pres - p1) / (p2 - p1);
                    flag = 1;
                }
            }
        }

        /* Fall through - if we still didn't find a suitable pressure layer, then set the
		output value to be the highest value of the coarse array */
        if (retv == 0) {
            retv = y_coarse[yc_count];
        }

        return (retv);
    }

    /* This returns the dew point temperature (in KELVIN)
	given the specific humidity (q) in g/g
	and the pressure (p) in millibars. */
    /**
     * Method declaration
     *
     *
     * @param sp_hum
     * @param pres
     *
     * @return
     *
     * @see
     */
    public double _dpt_from_qp(double sp_hum, double pres) {

        double ew;
        double tdew;

        ew = sp_hum * pres / 0.622;
        tdew = 1.0 / ((1 / 273.15) - (461.51 / 2.49e6) * Math.log(ew / 6.11));

        return (tdew);

    }

    /*
	Calclulate the wet bulb potential temperature
     */
    public double _wetblb(double ppp, double ttt, double tttd) {

        /* Some constants */
        double r = 0.28704;
        double cp = 1.005;
        double rv = 0.46213;
        double xkap = 0.286;

        /* Declare local variables */
        double p, t, td, esat, e, theta;
        /* Latent heat of evaporation */
        double rns, rs, pllv, pllvrs, dt;

        double dp = 160;
        double p0 = 1000.0;

        double thw = 0;

        /* Preserve input parameters */
        p = ppp;
        t = ttt;
        td = tttd;

        /* Calculate mixing ratios (rns, rs) and potential temperature (theta) */
        esat = _satvap(t, td);
        e = _satvap(td, td);
        rns = _wfun(esat, p);
        rs = _wfun(e, p);
        theta = t * Math.pow((p0 / p), xkap);

        /* Decrement p by 80mb. Calculate new t using constant theta.
		loop until rns<rs. This occurs only after p < pcl. */
        if ((ttt - tttd) > 0.01) /* Skip calculation if tttd~=ttt */ {
            while (rns > rs) {
                p = p - 80.0;
                if (p <= 0) {
                    p = 0.00001;
                }
                t = theta / Math.pow((p0 / p), xkap);
                rns = _wfun(_satvap(t, t), p);
            }

            /* Now zero in on pressure lifing condensation level until criterion in first "if" is met. */
            while (!((Math.abs(rns - rs) <= 5) | (Math.abs(dp) < 1))) // Ask TNC expr eval to 5e -5
            {
                dp = dp / 2.0;
                /* Reduce pressure step */
                if ((rns - rs) < 0) {
                    p = p + dp;
                }
                if ((rns - rs) >= 0) {
                    p = p - dp;
                }
                t = theta / Math.pow(p0 / p, xkap);
                rns = _wfun(_satvap(t, t), p);
            }
        }

        /* Go down moist adiabat to 1000MB or sfc press. loop until criterion is met.
		See the ceaseless wind for dt calculation (278,16) */
        dp = 20.0;

        do {
            if (p > (p0 - dp)) {
                dp = p0 - p;
            }
            p = p + 0.5 * dp;
            rs = _wfun(_satvap(t, t), p);
            pllv = _levap(t);
            pllvrs = pllv * rs;
            dt = (((r + pllvrs / t) * (dp / p)) / (cp + (pllv * pllvrs) / (rv * t * t))) * t;
            p = p + 0.5 * dp;
            t = t + dt;
        } while (Math.abs(p - p0) > 0.5);
        thw = t;
        return (thw);	// Return the wet bulb  potential temperature.
    }


    /* Some small functions called by the wet bulb
	potential temperature routines. I don't know
	what these are calculating, because the fortran
	source doesn't deign to say. Maybe it's obvious
	to a meteorologist. */
    public double _levap(double t) {
        return ((597.3 - 0.57 * (t - 273.17)) * 4.186);
    }

    public double _satvap(double td, double t) {
        double rv = 0.46213;
        return (6.11 * Math.exp(_levap(t) * (1.0 / 273.0 - 1.0 / td) / rv));
    }

    public double _wfun(double e, double p) {
        return (0.622 * e / (p - e));
    }

    /* Computes thte from pres, tmpc, dwpc. In the 
        calculation, mixr depends on pres and dwpc; tlcl
        depends on tmpc and dwpc. The following equation is
        used (deleted, can't be bothered.)

        Input parameters:
        pres -> pressure in millibars;
        tmpc -> temperature in celcius;
        dwpc -> dewpoint in celcius.

        output parameter:
        pr_thte -> equivalent potential temperature in K. */
    public double _pr_thte(double pres, double tmpc, double dwpc) {

        double pr_thte, rmix, tmpk, dwpk, e, thtam, tlcl, vapr, corr, eth;

        /* Check for missing values, except really we don't */
 /* Find mixing ratio; check for bad values */
 /* Function to compute mixr from dwpc and pres. The following equation is used:
                mixr = 0.62197 * (e/(pres - e)) * 1000.0;
                e = vapr * corr;
                corr = (1.001 + ((pres - 100.0)/900.0) * 0.0034)
                University of Wisconsin Green Sheet.
         */
 /* Calculate vapour pressure */
        vapr = 6.112 * Math.exp((17.67 * dwpc) / (dwpc + 243.5));

        /* Corr is a correction to the vapour pressure
                since the atmosphere is not an ideal gas. */
        corr = (1.001 + ((pres - 100.0) / 900) * 0.0034);
        e = corr * vapr;

        /* Test for unphysical case of large e at low pres. */
 /*
                if(e > (0.5 * pres)) {
                        pr_mixr = RMISSD;
                }
                else {
         */
 /* Calculate mixing ratio */
        rmix = 0.62197 * (e / (pres - e)) * 1000.0;

        //}
        /* Change degress celcius to kelvin */
        tmpk = (tmpc + C_TO_K);

        /* Calclulate theta - m (theta for moist air) */
        eth = (2.0 / 7.0) * (1.0 - (0.28 * 0.001 * rmix));

        thtam = tmpk * Math.pow((1000.0 / pres), e);

        /* Find the temperature at the LCL */
 /* Function to compute temperature at the lifted condensation level
                for a parcel of air given tmpc and dwpc. The following equation is
                used:tlcl = (1/(1/(dwpk - 56) + alog(tmpk/dwpk)/800)) + 56
         */
        dwpk = (dwpc + C_TO_K);
        tlcl = (1.0 / (1.0 / (dwpk - 56.0) + Math.log(tmpk / dwpk) / 800.0)) + 56.0;

        e = ((3.376 / tlcl) - 0.00254) * (rmix * (1.0 + 0.81 * 0.001 * rmix));
        pr_thte = thtam * Math.exp(e);
        return (pr_thte);
    }

    /* This function computes tmst from thte, pres, tguess.
        tmst is the parcel temperature at level pres on a 
        specified moist adiabat (thte.) The computation is
        an iterative Newton-Rapheson technique of the form
        x = x(guess) + (f(x) - f(x(guess)))/f'(x(guess));
        f' is approximated with finite differences.
        f' = (f(x(guess) + 1) - f(x(guess)))/1

        Convergence is not guaranteed for extreme input values.
        If the computation does not converge after 100 iterations,
        the missing data value will be returned.

        input parameters:
        thte -> Equivalent potential temperatur in K.
        pres -> pressure in millibars.
        tguess -> first guess temperature in K.

        ouput parameter:
        pr_tmst-> Parcel temperature in kelvin. */
    public double _pr_tmst(double thte, double pres, double tguess) {
        double pr_tmst, tg, epsi, tgnu, tgnup, tenu, tenup, cor;
        int i;
        pr_tmst = 0;
        /*
                if(_ermiss(thte)||_ermiss(pres)||_ermiss(tguess)||(thte<=0)||(pres<=0)||(tguess<0))
                {
                        pr_tmst = RMISSD;
                }
                else
                {}
         */

 /* Move tguess into another variable. */
        tg = tguess;

        /* If tguess is passed as 0.0 it is computed from an MIT scheme. */
        if (tg == 0.0) {
            tg = (thte - 0.5 * Math.pow(_max(thte - 270.0, 0.0), 1.05)) * Math.pow(pres / 1000.0, 0.2);
        }

        /* Set convergence and initial guess in degress C */
        epsi = 0.01;
        tgnu = (tg - C_TO_K);

        /* Set a limit of 100 iterations. Computer tenu, tenup, the
                        thte's at, one degree above the guess temperature. */
        for (i = 1; i <= 100; i++) {
            tgnup = tgnu + 1.0;
            tenu = _pr_thte(pres, tgnu, tgnu);
            tenup = _pr_thte(pres, tgnup, tgnup);

            /* Check that the thte's exist. */
 /*
                        if(_ermiss(tenu)||_ermiss(tenup)){
                                pr_tmst = RMISSD;
                        }
                        else
                        {}
             */
 /* Compute the correction, deltg; return on convergence */
            cor = (thte - tenu) / (tenup - tenu);
            tgnu = tgnu + cor;
            if ((cor < epsi) && (-cor < epsi)) {
                pr_tmst = (tgnu + C_TO_K);
            }
        }
        /* Failed to converge - return missing */
 /*
                if(i==100) {
                        pr_tmst = RMISSD;
                }
         */
        return (pr_tmst);
    }

    /*
        public double _ermiss(double x) {
                return (Math.abs(x - RMISSD) < 0.1);
        }
     */
 /*
	 * Calculates the conductivity of the soil
	 * Cosby curves and coefficients.
     */
    /**
     * Method declaration
     *
     *
     * @param thmax
     * @param thv
     * @param psis
     * @param cosbyb
     *
     * @return
     *
     * @see
     */
    public double _psgcal(double thmax, double thv, double psis,
            double cosbyb) {
        double perwmax, perw2g, rlogpsig, rnpsig;

        /*
		 * Field capacity (75% of thmax) used instead of thmax.
		 * Lower value felt to fit local measurements better than
		 * the fit with tabulated values.
         */
 /* Convert ground water contents to percentages. */
        perwmax = (thmax * 100 * 0.75);
        perw2g = thv * 100;
        rlogpsig = log10(psis) + cosbyb * log10(perwmax)
                - cosbyb * log10(perw2g);

        /* psig is positive in this program. */
        rnpsig = (Math.pow(10, rlogpsig));

        /* Convert cm to bars and return psig */
        return ((-rnpsig / 1020));
    }

    /**
     * Method declaration
     *
     *
     * @param x
     *
     * @return
     *
     * @see
     */
    public double log10(double x) {
        return logBaseAofX(10, x);
    }

    /**
     * Method declaration
     *
     *
     * @param a
     * @param x
     *
     * @return
     *
     * @see
     */
    public double logBaseAofX(int a, double x) {

        // log base a (x) = log base b (x) / log base b (a)
        // = ln(x) / ln(a)
        return Math.log(x) / Math.log(a);
    }

    /**
     * Method declaration
     *
     *
     * @param top
     * @param bot
     * @param index
     *
     * @return
     *
     * @see
     */
    public double _rlinear(double top, double bot, int index) {
        return ((bot + (top - bot) * index / 5));
    }

    /*
	 * If below top of the mixing layer potential temperature
	 * that at 50m. Otherwise, average as usual to create fine scale.
     */
    /**
     * Method declaration
     *
     *
     * @param u_fine
     * @param v_fine
     * @param t_fine
     * @param q_fine
     * @param ud
     * @param vd
     * @param td
     * @param qd
     * @param hgt
     * @param aptemp
     *
     * @see
     */
    public void _fine(double[] u_fine, double[] v_fine, double[] t_fine,
            double[] q_fine, double[] ud, double[] vd, double[] td,
            double[] qd, double hgt, double aptemp) {
        int i, incr, j, height;

        height = 50;
        j = 1;
        u_fine[1] = ud[1];
        v_fine[1] = vd[1];
        t_fine[1] = td[1];
        q_fine[1] = qd[1];

        for (i = 1; i <= 9; i++) {
            for (incr = 1; incr <= 5; incr++) {
                j = j + 1;
                u_fine[j] = _rlinear(ud[i + 1], ud[i], incr);
                v_fine[j] = _rlinear(vd[i + 1], vd[i], incr);
                q_fine[j] = _rlinear(qd[i + 1], qd[i], incr);

                if (height <= hgt) {
                    t_fine[j] = aptemp;
                } else {
                    t_fine[j] = _rlinear(td[i + 1], td[i], incr);
                }

                height += 50;
            }
        }
    }

    // From spline.h
    public final int NR_END = 1;

    // #define FREE_ARG char*

    /* Cubic spline interpolation functions, from "Numerical Recipes in C" */
    /**
     * Method declaration
     *
     *
     * @param x
     * @param y
     * @param n
     * @param yp1
     * @param ypn
     * @param y2
     *
     * @see
     */
    public void spline(double x[], double y[], int n, double yp1, double ypn,
            double y2[]) {
        int i, k;
        double p, qn, sig, un;			 // , *u;

        // u=vector(1,n-1);
        double[] u
                = new double[(n - 1) - 1 + 1 + NR_END];	 // tcs see vector routine

        if (yp1 > 0.99e30) {
            y2[1] = (u[1] = 0.0);
        } else {
            y2[1] = -0.5;
            u[1] = ((3.0 / (x[2] - x[1]))
                    * ((y[2] - y[1]) / (x[2] - x[1]) - yp1));
        }

        for (i = 2; i <= n - 1; i++) {
            sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
            p = (sig * y2[i - 1] + 2.0);
            y2[i] = ((sig - 1.0) / p);
            u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
                    - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
            u[i] = ((6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1])
                    / p);
        }

        if (ypn > 0.99e30) {
            qn = un = 0.0;
        } else {
            qn = 0.5;
            un = ((3.0 / (x[n] - x[n - 1]))
                    * (ypn - (y[n] - y[n - 1]) / (x[n] - x[n - 1])));
        }

        y2[n] = ((un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0));

        for (k = n - 1; k >= 1; k--) {
            y2[k] = y2[k] * y2[k + 1] + u[k];
        }

        // free_vector(u,1,n-1);
        // u should garbage collect eventually
        // System.gc();
        // System.gc();
        // System.gc();
    }

    // void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
    // tcs passes in a ptr to a double, then seems to be just a scalar at the end, so
    // I changed it to return a double
    /**
     * Method declaration
     *
     *
     * @param xa
     * @param ya
     * @param y2a
     * @param n
     * @param x
     *
     * @return
     *
     * @throws Exception
     *
     * @see
     */
    public double splint(double xa[], double ya[], double y2a[], int n,
            double x) throws Exception {

        // void nrerror(char error_text[]);
        // tcs wacky error catcher

        /*
		 * TCS Attempt to work around JIT error - what a pain
		 * Here's the original code...
		 * 
		 * int klo, khi, k;
		 * double h, b, a;
		 * 
		 * klo=1;
		 * khi=n;
		 * 
		 * while(khi-klo>1)
		 * {
		 * k=(khi+klo)>>1;
		 * if(xa[k]>x)
		 * {
		 * khi=k;
		 * }
		 * else
		 * {
		 * klo=k;
		 * }
		 * }
		 * 
		 * h=xa[khi]-xa[klo];
		 * if(h==0)
		 * {
		 * //nrerror("Bad xa input to routine splint");
		 * Exception t = new Exception("Bad xa input to routine splint");
		 * throw t;
		 * }
		 * a=(xa[khi]-x)/h;
		 * b=(x-xa[klo])/h;
		 * return (a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0);
         */
        int klo, khi, k;
        double h, b, a;

        klo = 1;
        khi = n;

        while ((khi - klo) > 1) {
            k = khi + klo;	  // TCS I don't think it liked k = (khi+klo) >> 1;
            k = k >> 1;

            if (xa[k] > x) {
                khi = k;
            } else {
                klo = k;
            }
        }

        h = xa[khi] - xa[klo];

        if (h == 0) {
            Exception t = new Exception("Bad xa input to routine splint");

            throw t;
        }

        a = (xa[khi] - x) / h;
        b = (x - xa[klo]) / h;

        // TCS changed a*a*a to Math.pow(a,3) and b*b*b to Math.pow(b,3)
        return a * ya[klo] + b * ya[khi]
                + ((Math.pow(a, 3) - a) * y2a[klo] + (Math.pow(b, 3) - b) * y2a[khi])
                * (h * h) / 6.0;
    }

    /*
	 * tcs error routine - we have try/catch
	 * void nrerror(char error_text[])
	 * {
	 * System.out.println("Numerical recipes run-time error...\n");
	 * System.out.println("%s\n",error_text);
	 * System.out.println("Exiting...\n");
	 * throw new Exception();
	 * }
     */
    // tcs Some wacky double array. I'll try just new double[]

    /*
	 * double *vector(long nl, long nh)
	 * {
	 * double *v;
	 * v=(double *)malloc((size_t)((nh-nl+1+NR_END)*sizeof));
	 * if(!v)
	 * {
	 * nrerror("Allocation failure in vector()");
	 * }
	 * return(v-nl+NR_END);
	 * }
     */

 /*
	 * void free_vector(double *v, long nl, long nh)
	 * {
	 * free((FREE_ARG)(v+nl-NR_END));
	 * }
     */
}



/*--- formatting done in "Sun Java Convention" style on 02-15-2000 ---*/
