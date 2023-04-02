//
// Created by wcy on 2023/2/27.
//
#include <math.h>
#include <malloc.h>
#include <stdio.h>

#include "real_fft.h"

static void radb2(int ido, int l1, const real *cc, real *ch, const real *wa1) {
    /* System generated locals */
    int cc_offset, ch_offset;

    /* Local variables */
    int i, k, ic;
    real ti2, tr2;
    int idp2;


#define cc_ref(a_1, a_2, a_3) cc[((a_3)*2 + (a_2))*ido + a_1]
#define ch_ref(a_1, a_2, a_3) ch[((a_3)*l1 + (a_2))*ido + a_1]

    /* Parameter adjustments */
    ch_offset = 1 + ido * (1 + l1);
    ch -= ch_offset;
    cc_offset = 1 + ido * 3;
    cc -= cc_offset;
    --wa1;

    /* Function Body */
    for (k = 1; k <= l1; ++k) {
        ch_ref(1, k, 1) = cc_ref(1, 1, k) + cc_ref(ido, 2, k);
        ch_ref(1, k, 2) = cc_ref(1, 1, k) - cc_ref(ido, 2, k);
    }
    if (ido < 2) return;
    else if (ido != 2) {
        idp2 = ido + 2;
        for (k = 1; k <= l1; ++k) {
            for (i = 3; i <= ido; i += 2) {
                ic = idp2 - i;
                ch_ref(i - 1, k, 1) = cc_ref(i - 1, 1, k) + cc_ref(ic - 1, 2,
                                                                   k);
                tr2 = cc_ref(i - 1, 1, k) - cc_ref(ic - 1, 2, k);
                ch_ref(i, k, 1) = cc_ref(i, 1, k) - cc_ref(ic, 2, k);
                ti2 = cc_ref(i, 1, k) + cc_ref(ic, 2, k);
                ch_ref(i - 1, k, 2) = wa1[i - 2] * tr2 - wa1[i - 1] * ti2;
                ch_ref(i, k, 2) = wa1[i - 2] * ti2 + wa1[i - 1] * tr2;
            }
        }
        if (ido % 2 == 1) return;
    }
    for (k = 1; k <= l1; ++k) {
        ch_ref(ido, k, 1) = cc_ref(ido, 1, k) + cc_ref(ido, 1, k);
        ch_ref(ido, k, 2) = -(cc_ref(1, 2, k) + cc_ref(1, 2, k));
    }
} /* radb2 */

#undef ch_ref
#undef cc_ref


static void radb3(int ido, int l1, const real *cc, real *ch,
                  const real *wa1, const real *wa2) {
    /* Initialized data */

    static const real taur = -.5f;
    static const real taui = .866025403784439f;

    /* System generated locals */
    int cc_offset, ch_offset;

    /* Local variables */
    int i, k, ic;
    real ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
    int idp2;


#define cc_ref(a_1, a_2, a_3) cc[((a_3)*3 + (a_2))*ido + a_1]
#define ch_ref(a_1, a_2, a_3) ch[((a_3)*l1 + (a_2))*ido + a_1]

    /* Parameter adjustments */
    ch_offset = 1 + ido * (1 + l1);
    ch -= ch_offset;
    cc_offset = 1 + (ido << 2);
    cc -= cc_offset;
    --wa1;
    --wa2;

    /* Function Body */
    for (k = 1; k <= l1; ++k) {
        tr2 = cc_ref(ido, 2, k) + cc_ref(ido, 2, k);
        cr2 = cc_ref(1, 1, k) + taur * tr2;
        ch_ref(1, k, 1) = cc_ref(1, 1, k) + tr2;
        ci3 = taui * (cc_ref(1, 3, k) + cc_ref(1, 3, k));
        ch_ref(1, k, 2) = cr2 - ci3;
        ch_ref(1, k, 3) = cr2 + ci3;
    }
    if (ido == 1) {
        return;
    }
    idp2 = ido + 2;
    for (k = 1; k <= l1; ++k) {
        for (i = 3; i <= ido; i += 2) {
            ic = idp2 - i;
            tr2 = cc_ref(i - 1, 3, k) + cc_ref(ic - 1, 2, k);
            cr2 = cc_ref(i - 1, 1, k) + taur * tr2;
            ch_ref(i - 1, k, 1) = cc_ref(i - 1, 1, k) + tr2;
            ti2 = cc_ref(i, 3, k) - cc_ref(ic, 2, k);
            ci2 = cc_ref(i, 1, k) + taur * ti2;
            ch_ref(i, k, 1) = cc_ref(i, 1, k) + ti2;
            cr3 = taui * (cc_ref(i - 1, 3, k) - cc_ref(ic - 1, 2, k));
            ci3 = taui * (cc_ref(i, 3, k) + cc_ref(ic, 2, k));
            dr2 = cr2 - ci3;
            dr3 = cr2 + ci3;
            di2 = ci2 + cr3;
            di3 = ci2 - cr3;
            ch_ref(i - 1, k, 2) = wa1[i - 2] * dr2 - wa1[i - 1] * di2;
            ch_ref(i, k, 2) = wa1[i - 2] * di2 + wa1[i - 1] * dr2;
            ch_ref(i - 1, k, 3) = wa2[i - 2] * dr3 - wa2[i - 1] * di3;
            ch_ref(i, k, 3) = wa2[i - 2] * di3 + wa2[i - 1] * dr3;
        }
    }
} /* radb3 */

#undef ch_ref
#undef cc_ref


static void radb4(int ido, int l1, const real *cc, real *ch,
                  const real *wa1, const real *wa2, const real *wa3) {
    /* Initialized data */

    static const real sqrt2 = 1.414213562373095f;

    /* System generated locals */
    int cc_offset, ch_offset;

    /* Local variables */
    int i, k, ic;
    real ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
    int idp2;


#define cc_ref(a_1, a_2, a_3) cc[((a_3)*4 + (a_2))*ido + a_1]
#define ch_ref(a_1, a_2, a_3) ch[((a_3)*l1 + (a_2))*ido + a_1]

    /* Parameter adjustments */
    ch_offset = 1 + ido * (1 + l1);
    ch -= ch_offset;
    cc_offset = 1 + ido * 5;
    cc -= cc_offset;
    --wa1;
    --wa2;
    --wa3;

    /* Function Body */
    for (k = 1; k <= l1; ++k) {
        tr1 = cc_ref(1, 1, k) - cc_ref(ido, 4, k);
        tr2 = cc_ref(1, 1, k) + cc_ref(ido, 4, k);
        tr3 = cc_ref(ido, 2, k) + cc_ref(ido, 2, k);
        tr4 = cc_ref(1, 3, k) + cc_ref(1, 3, k);
        ch_ref(1, k, 1) = tr2 + tr3;
        ch_ref(1, k, 2) = tr1 - tr4;
        ch_ref(1, k, 3) = tr2 - tr3;
        ch_ref(1, k, 4) = tr1 + tr4;
    }
    if (ido < 2) return;
    if (ido != 2) {
        idp2 = ido + 2;
        for (k = 1; k <= l1; ++k) {
            for (i = 3; i <= ido; i += 2) {
                ic = idp2 - i;
                ti1 = cc_ref(i, 1, k) + cc_ref(ic, 4, k);
                ti2 = cc_ref(i, 1, k) - cc_ref(ic, 4, k);
                ti3 = cc_ref(i, 3, k) - cc_ref(ic, 2, k);
                tr4 = cc_ref(i, 3, k) + cc_ref(ic, 2, k);
                tr1 = cc_ref(i - 1, 1, k) - cc_ref(ic - 1, 4, k);
                tr2 = cc_ref(i - 1, 1, k) + cc_ref(ic - 1, 4, k);
                ti4 = cc_ref(i - 1, 3, k) - cc_ref(ic - 1, 2, k);
                tr3 = cc_ref(i - 1, 3, k) + cc_ref(ic - 1, 2, k);
                ch_ref(i - 1, k, 1) = tr2 + tr3;
                cr3 = tr2 - tr3;
                ch_ref(i, k, 1) = ti2 + ti3;
                ci3 = ti2 - ti3;
                cr2 = tr1 - tr4;
                cr4 = tr1 + tr4;
                ci2 = ti1 + ti4;
                ci4 = ti1 - ti4;
                ch_ref(i - 1, k, 2) = wa1[i - 2] * cr2 - wa1[i - 1] * ci2;
                ch_ref(i, k, 2) = wa1[i - 2] * ci2 + wa1[i - 1] * cr2;
                ch_ref(i - 1, k, 3) = wa2[i - 2] * cr3 - wa2[i - 1] * ci3;
                ch_ref(i, k, 3) = wa2[i - 2] * ci3 + wa2[i - 1] * cr3;
                ch_ref(i - 1, k, 4) = wa3[i - 2] * cr4 - wa3[i - 1] * ci4;
                ch_ref(i, k, 4) = wa3[i - 2] * ci4 + wa3[i - 1] * cr4;
            }
        }
        if (ido % 2 == 1) return;
    }
    for (k = 1; k <= l1; ++k) {
        ti1 = cc_ref(1, 2, k) + cc_ref(1, 4, k);
        ti2 = cc_ref(1, 4, k) - cc_ref(1, 2, k);
        tr1 = cc_ref(ido, 1, k) - cc_ref(ido, 3, k);
        tr2 = cc_ref(ido, 1, k) + cc_ref(ido, 3, k);
        ch_ref(ido, k, 1) = tr2 + tr2;
        ch_ref(ido, k, 2) = sqrt2 * (tr1 - ti1);
        ch_ref(ido, k, 3) = ti2 + ti2;
        ch_ref(ido, k, 4) = -sqrt2 * (tr1 + ti1);
    }
} /* radb4 */

#undef ch_ref
#undef cc_ref


static void radb5(int ido, int l1, const real *cc, real *ch,
                  const real *wa1, const real *wa2, const real *wa3, const real *wa4) {
    /* Initialized data */

    static const real tr11 = .309016994374947f;
    static const real ti11 = .951056516295154f;
    static const real tr12 = -.809016994374947f;
    static const real ti12 = .587785252292473f;

    /* System generated locals */
    int cc_offset, ch_offset;

    /* Local variables */
    int i, k, ic;
    real ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3,
            ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;
    int idp2;


#define cc_ref(a_1, a_2, a_3) cc[((a_3)*5 + (a_2))*ido + a_1]
#define ch_ref(a_1, a_2, a_3) ch[((a_3)*l1 + (a_2))*ido + a_1]

    /* Parameter adjustments */
    ch_offset = 1 + ido * (1 + l1);
    ch -= ch_offset;
    cc_offset = 1 + ido * 6;
    cc -= cc_offset;
    --wa1;
    --wa2;
    --wa3;
    --wa4;

    /* Function Body */
    for (k = 1; k <= l1; ++k) {
        ti5 = cc_ref(1, 3, k) + cc_ref(1, 3, k);
        ti4 = cc_ref(1, 5, k) + cc_ref(1, 5, k);
        tr2 = cc_ref(ido, 2, k) + cc_ref(ido, 2, k);
        tr3 = cc_ref(ido, 4, k) + cc_ref(ido, 4, k);
        ch_ref(1, k, 1) = cc_ref(1, 1, k) + tr2 + tr3;
        cr2 = cc_ref(1, 1, k) + tr11 * tr2 + tr12 * tr3;
        cr3 = cc_ref(1, 1, k) + tr12 * tr2 + tr11 * tr3;
        ci5 = ti11 * ti5 + ti12 * ti4;
        ci4 = ti12 * ti5 - ti11 * ti4;
        ch_ref(1, k, 2) = cr2 - ci5;
        ch_ref(1, k, 3) = cr3 - ci4;
        ch_ref(1, k, 4) = cr3 + ci4;
        ch_ref(1, k, 5) = cr2 + ci5;
    }
    if (ido == 1) {
        return;
    }
    idp2 = ido + 2;
    for (k = 1; k <= l1; ++k) {
        for (i = 3; i <= ido; i += 2) {
            ic = idp2 - i;
            ti5 = cc_ref(i, 3, k) + cc_ref(ic, 2, k);
            ti2 = cc_ref(i, 3, k) - cc_ref(ic, 2, k);
            ti4 = cc_ref(i, 5, k) + cc_ref(ic, 4, k);
            ti3 = cc_ref(i, 5, k) - cc_ref(ic, 4, k);
            tr5 = cc_ref(i - 1, 3, k) - cc_ref(ic - 1, 2, k);
            tr2 = cc_ref(i - 1, 3, k) + cc_ref(ic - 1, 2, k);
            tr4 = cc_ref(i - 1, 5, k) - cc_ref(ic - 1, 4, k);
            tr3 = cc_ref(i - 1, 5, k) + cc_ref(ic - 1, 4, k);
            ch_ref(i - 1, k, 1) = cc_ref(i - 1, 1, k) + tr2 + tr3;
            ch_ref(i, k, 1) = cc_ref(i, 1, k) + ti2 + ti3;
            cr2 = cc_ref(i - 1, 1, k) + tr11 * tr2 + tr12 * tr3;
            ci2 = cc_ref(i, 1, k) + tr11 * ti2 + tr12 * ti3;
            cr3 = cc_ref(i - 1, 1, k) + tr12 * tr2 + tr11 * tr3;
            ci3 = cc_ref(i, 1, k) + tr12 * ti2 + tr11 * ti3;
            cr5 = ti11 * tr5 + ti12 * tr4;
            ci5 = ti11 * ti5 + ti12 * ti4;
            cr4 = ti12 * tr5 - ti11 * tr4;
            ci4 = ti12 * ti5 - ti11 * ti4;
            dr3 = cr3 - ci4;
            dr4 = cr3 + ci4;
            di3 = ci3 + cr4;
            di4 = ci3 - cr4;
            dr5 = cr2 + ci5;
            dr2 = cr2 - ci5;
            di5 = ci2 - cr5;
            di2 = ci2 + cr5;
            ch_ref(i - 1, k, 2) = wa1[i - 2] * dr2 - wa1[i - 1] * di2;
            ch_ref(i, k, 2) = wa1[i - 2] * di2 + wa1[i - 1] * dr2;
            ch_ref(i - 1, k, 3) = wa2[i - 2] * dr3 - wa2[i - 1] * di3;
            ch_ref(i, k, 3) = wa2[i - 2] * di3 + wa2[i - 1] * dr3;
            ch_ref(i - 1, k, 4) = wa3[i - 2] * dr4 - wa3[i - 1] * di4;
            ch_ref(i, k, 4) = wa3[i - 2] * di4 + wa3[i - 1] * dr4;
            ch_ref(i - 1, k, 5) = wa4[i - 2] * dr5 - wa4[i - 1] * di5;
            ch_ref(i, k, 5) = wa4[i - 2] * di5 + wa4[i - 1] * dr5;
        }
    }
} /* radb5 */

#undef ch_ref
#undef cc_ref


static void radbg(int ido, int ip, int l1, int idl1,
                  const real *cc, real *c1, real *c2, real *ch, real *ch2, const real *wa) {
    /* System generated locals */
    int ch_offset, cc_offset,
            c1_offset, c2_offset, ch2_offset;

    /* Local variables */
    int i, j, k, l, j2, ic, jc, lc, ik, is;
    real dc2, ai1, ai2, ar1, ar2, ds2;
    int nbd;
    real dcp, arg, dsp, ar1h, ar2h;
    int idp2, ipp2, idij, ipph;


#define c1_ref(a_1, a_2, a_3) c1[((a_3)*l1 + (a_2))*ido + a_1]
#define c2_ref(a_1, a_2) c2[(a_2)*idl1 + a_1]
#define cc_ref(a_1, a_2, a_3) cc[((a_3)*ip + (a_2))*ido + a_1]
#define ch_ref(a_1, a_2, a_3) ch[((a_3)*l1 + (a_2))*ido + a_1]
#define ch2_ref(a_1, a_2) ch2[(a_2)*idl1 + a_1]

    ch_offset = 1 + ido * (1 + l1);
    ch -= ch_offset;
    c1_offset = 1 + ido * (1 + l1);
    c1 -= c1_offset;
    cc_offset = 1 + ido * (1 + ip);
    cc -= cc_offset;
    ch2_offset = 1 + idl1;
    ch2 -= ch2_offset;
    c2_offset = 1 + idl1;
    c2 -= c2_offset;
    --wa;

    arg = (2 * M_PI) / (real) (ip);
    dcp = real_cos(arg);
    dsp = real_sin(arg);
    idp2 = ido + 2;
    nbd = (ido - 1) / 2;
    ipp2 = ip + 2;
    ipph = (ip + 1) / 2;
    if (ido >= l1) {
        for (k = 1; k <= l1; ++k) {
            for (i = 1; i <= ido; ++i) {
                ch_ref(i, k, 1) = cc_ref(i, 1, k);
            }
        }
    } else {
        for (i = 1; i <= ido; ++i) {
            for (k = 1; k <= l1; ++k) {
                ch_ref(i, k, 1) = cc_ref(i, 1, k);
            }
        }
    }
    for (j = 2; j <= ipph; ++j) {
        jc = ipp2 - j;
        j2 = j + j;
        for (k = 1; k <= l1; ++k) {
            ch_ref(1, k, j) = cc_ref(ido, j2 - 2, k) + cc_ref(ido, j2 - 2, k);
            ch_ref(1, k, jc) = cc_ref(1, j2 - 1, k) + cc_ref(1, j2 - 1, k);
        }
    }
    if (ido != 1) {
        if (nbd >= l1) {
            for (j = 2; j <= ipph; ++j) {
                jc = ipp2 - j;
                for (k = 1; k <= l1; ++k) {
                    for (i = 3; i <= ido; i += 2) {
                        ic = idp2 - i;
                        ch_ref(i - 1, k, j) = cc_ref(i - 1, (j << 1) - 1, k) + cc_ref(ic - 1, (j << 1) - 2, k);
                        ch_ref(i - 1, k, jc) = cc_ref(i - 1, (j << 1) - 1, k) - cc_ref(ic - 1, (j << 1) - 2, k);
                        ch_ref(i, k, j) = cc_ref(i, (j << 1) - 1, k) - cc_ref(ic, (j << 1) - 2, k);
                        ch_ref(i, k, jc) = cc_ref(i, (j << 1) - 1, k) + cc_ref(ic, (j << 1) - 2, k);
                    }
                }
            }
        } else {
            for (j = 2; j <= ipph; ++j) {
                jc = ipp2 - j;
                for (i = 3; i <= ido; i += 2) {
                    ic = idp2 - i;
                    for (k = 1; k <= l1; ++k) {
                        ch_ref(i - 1, k, j) = cc_ref(i - 1, (j << 1) - 1, k) + cc_ref(ic - 1, (j << 1) - 2, k);
                        ch_ref(i - 1, k, jc) = cc_ref(i - 1, (j << 1) - 1, k) - cc_ref(ic - 1, (j << 1) - 2, k);
                        ch_ref(i, k, j) = cc_ref(i, (j << 1) - 1, k) - cc_ref(ic, (j << 1) - 2, k);
                        ch_ref(i, k, jc) = cc_ref(i, (j << 1) - 1, k) + cc_ref(ic, (j << 1) - 2, k);
                    }
                }
            }
        }
    }
    ar1 = 1.f;
    ai1 = 0.f;
    for (l = 2; l <= ipph; ++l) {
        lc = ipp2 - l;
        ar1h = dcp * ar1 - dsp * ai1;
        ai1 = dcp * ai1 + dsp * ar1;
        ar1 = ar1h;
        for (ik = 1; ik <= idl1; ++ik) {
            c2_ref(ik, l) = ch2_ref(ik, 1) + ar1 * ch2_ref(ik, 2);
            c2_ref(ik, lc) = ai1 * ch2_ref(ik, ip);
        }
        dc2 = ar1;
        ds2 = ai1;
        ar2 = ar1;
        ai2 = ai1;
        for (j = 3; j <= ipph; ++j) {
            jc = ipp2 - j;
            ar2h = dc2 * ar2 - ds2 * ai2;
            ai2 = dc2 * ai2 + ds2 * ar2;
            ar2 = ar2h;
            for (ik = 1; ik <= idl1; ++ik) {
                c2_ref(ik, l) = c2_ref(ik, l) + ar2 * ch2_ref(ik, j);
                c2_ref(ik, lc) = c2_ref(ik, lc) + ai2 * ch2_ref(ik, jc);
            }
        }
    }
    for (j = 2; j <= ipph; ++j) {
        for (ik = 1; ik <= idl1; ++ik) {
            ch2_ref(ik, 1) = ch2_ref(ik, 1) + ch2_ref(ik, j);
        }
    }
    for (j = 2; j <= ipph; ++j) {
        jc = ipp2 - j;
        for (k = 1; k <= l1; ++k) {
            ch_ref(1, k, j) = c1_ref(1, k, j) - c1_ref(1, k, jc);
            ch_ref(1, k, jc) = c1_ref(1, k, j) + c1_ref(1, k, jc);
        }
    }
    if (ido != 1) {
        if (nbd >= l1) {
            for (j = 2; j <= ipph; ++j) {
                jc = ipp2 - j;
                for (k = 1; k <= l1; ++k) {
                    for (i = 3; i <= ido; i += 2) {
                        ch_ref(i - 1, k, j) = c1_ref(i - 1, k, j) - c1_ref(i, k, jc);
                        ch_ref(i - 1, k, jc) = c1_ref(i - 1, k, j) + c1_ref(i, k, jc);
                        ch_ref(i, k, j) = c1_ref(i, k, j) + c1_ref(i - 1, k, jc);
                        ch_ref(i, k, jc) = c1_ref(i, k, j) - c1_ref(i - 1, k, jc);
                    }
                }
            }
        } else {
            for (j = 2; j <= ipph; ++j) {
                jc = ipp2 - j;
                for (i = 3; i <= ido; i += 2) {
                    for (k = 1; k <= l1; ++k) {
                        ch_ref(i - 1, k, j) = c1_ref(i - 1, k, j) - c1_ref(i, k, jc);
                        ch_ref(i - 1, k, jc) = c1_ref(i - 1, k, j) + c1_ref(i, k, jc);
                        ch_ref(i, k, j) = c1_ref(i, k, j) + c1_ref(i - 1, k, jc);
                        ch_ref(i, k, jc) = c1_ref(i, k, j) - c1_ref(i - 1, k, jc);
                    }
                }
            }
        }
    }
    if (ido == 1) {
        return;
    }
    for (ik = 1; ik <= idl1; ++ik) {
        c2_ref(ik, 1) = ch2_ref(ik, 1);
    }
    for (j = 2; j <= ip; ++j) {
        for (k = 1; k <= l1; ++k) {
            c1_ref(1, k, j) = ch_ref(1, k, j);
        }
    }
    if (nbd <= l1) {
        is = -(ido);
        for (j = 2; j <= ip; ++j) {
            is += ido;
            idij = is;
            for (i = 3; i <= ido; i += 2) {
                idij += 2;
                for (k = 1; k <= l1; ++k) {
                    c1_ref(i - 1, k, j) = wa[idij - 1] * ch_ref(i - 1, k, j) - wa[idij] * ch_ref(i, k, j);
                    c1_ref(i, k, j) = wa[idij - 1] * ch_ref(i, k, j) + wa[idij] * ch_ref(i - 1, k, j);
                }
            }
        }
    } else {
        is = -(ido);
        for (j = 2; j <= ip; ++j) {
            is += ido;
            for (k = 1; k <= l1; ++k) {
                idij = is;
                for (i = 3; i <= ido; i += 2) {
                    idij += 2;
                    c1_ref(i - 1, k, j) = wa[idij - 1] * ch_ref(i - 1, k, j) - wa[idij] * ch_ref(i, k, j);
                    c1_ref(i, k, j) = wa[idij - 1] * ch_ref(i, k, j) + wa[idij] * ch_ref(i - 1, k, j);
                }
            }
        }
    }
}

#undef ch2_ref
#undef ch_ref
#undef cc_ref
#undef c2_ref
#undef c1_ref


static void radf2(int ido, int l1, const real *cc, real *ch,
                  const real *wa1) {
    int ch_offset, cc_offset;

    int i, k, ic;
    real ti2, tr2;
    int idp2;


#define cc_ref(a_1, a_2, a_3) cc[((a_3)*l1 + (a_2))*ido + a_1]
#define ch_ref(a_1, a_2, a_3) ch[((a_3)*2 + (a_2))*ido + a_1]

    ch_offset = 1 + ido * 3;
    ch -= ch_offset;
    cc_offset = 1 + ido * (1 + l1);
    cc -= cc_offset;
    --wa1;

    for (k = 1; k <= l1; ++k) {
        ch_ref(1, 1, k) = cc_ref(1, k, 1) + cc_ref(1, k, 2);
        ch_ref(ido, 2, k) = cc_ref(1, k, 1) - cc_ref(1, k, 2);
    }
    if (ido < 2) return;
    if (ido != 2) {
        idp2 = ido + 2;
        for (k = 1; k <= l1; ++k) {
            for (i = 3; i <= ido; i += 2) {
                ic = idp2 - i;
                tr2 = wa1[i - 2] * cc_ref(i - 1, k, 2) + wa1[i - 1] * cc_ref(i, k, 2);
                ti2 = wa1[i - 2] * cc_ref(i, k, 2) - wa1[i - 1] * cc_ref(i - 1, k, 2);
                ch_ref(i, 1, k) = cc_ref(i, k, 1) + ti2;
                ch_ref(ic, 2, k) = ti2 - cc_ref(i, k, 1);
                ch_ref(i - 1, 1, k) = cc_ref(i - 1, k, 1) + tr2;
                ch_ref(ic - 1, 2, k) = cc_ref(i - 1, k, 1) - tr2;
            }
        }
        if (ido % 2 == 1) {
            return;
        }
    }
    for (k = 1; k <= l1; ++k) {
        ch_ref(1, 2, k) = -cc_ref(ido, k, 2);
        ch_ref(ido, 1, k) = cc_ref(ido, k, 1);
    }
}

#undef ch_ref
#undef cc_ref


static void radf3(int ido, int l1, const real *cc, real *ch,
                  const real *wa1, const real *wa2) {
    static const real taur = -.5f;
    static const real taui = .866025403784439f;

    int ch_offset, cc_offset;

    int i, k, ic;
    real ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3;
    int idp2;


#define cc_ref(a_1, a_2, a_3) cc[((a_3)*l1 + (a_2))*ido + a_1]
#define ch_ref(a_1, a_2, a_3) ch[((a_3)*3 + (a_2))*ido + a_1]

    ch_offset = 1 + (ido << 2);
    ch -= ch_offset;
    cc_offset = 1 + ido * (1 + l1);
    cc -= cc_offset;
    --wa1;
    --wa2;

    for (k = 1; k <= l1; ++k) {
        cr2 = cc_ref(1, k, 2) + cc_ref(1, k, 3);
        ch_ref(1, 1, k) = cc_ref(1, k, 1) + cr2;
        ch_ref(1, 3, k) = taui * (cc_ref(1, k, 3) - cc_ref(1, k, 2));
        ch_ref(ido, 2, k) = cc_ref(1, k, 1) + taur * cr2;
    }
    if (ido == 1) {
        return;
    }
    idp2 = ido + 2;
    for (k = 1; k <= l1; ++k) {
        for (i = 3; i <= ido; i += 2) {
            ic = idp2 - i;
            dr2 = wa1[i - 2] * cc_ref(i - 1, k, 2) + wa1[i - 1] * cc_ref(i, k, 2);
            di2 = wa1[i - 2] * cc_ref(i, k, 2) - wa1[i - 1] * cc_ref(i - 1, k, 2);
            dr3 = wa2[i - 2] * cc_ref(i - 1, k, 3) + wa2[i - 1] * cc_ref(i, k, 3);
            di3 = wa2[i - 2] * cc_ref(i, k, 3) - wa2[i - 1] * cc_ref(i - 1, k, 3);
            cr2 = dr2 + dr3;
            ci2 = di2 + di3;
            ch_ref(i - 1, 1, k) = cc_ref(i - 1, k, 1) + cr2;
            ch_ref(i, 1, k) = cc_ref(i, k, 1) + ci2;
            tr2 = cc_ref(i - 1, k, 1) + taur * cr2;
            ti2 = cc_ref(i, k, 1) + taur * ci2;
            tr3 = taui * (di2 - di3);
            ti3 = taui * (dr3 - dr2);
            ch_ref(i - 1, 3, k) = tr2 + tr3;
            ch_ref(ic - 1, 2, k) = tr2 - tr3;
            ch_ref(i, 3, k) = ti2 + ti3;
            ch_ref(ic, 2, k) = ti3 - ti2;
        }
    }
}

#undef ch_ref
#undef cc_ref


static void radf4(int ido, int l1, const real *cc, real *ch,
                  const real *wa1, const real *wa2, const real *wa3) {
    static const real hsqt2 = .7071067811865475f;

    int cc_offset, ch_offset;

    int i, k, ic;
    real ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
    int idp2;


#define cc_ref(a_1, a_2, a_3) cc[((a_3)*l1 + (a_2))*ido + a_1]
#define ch_ref(a_1, a_2, a_3) ch[((a_3)*4 + (a_2))*ido + a_1]

    ch_offset = 1 + ido * 5;
    ch -= ch_offset;
    cc_offset = 1 + ido * (1 + l1);
    cc -= cc_offset;
    --wa1;
    --wa2;
    --wa3;

    for (k = 1; k <= l1; ++k) {
        tr1 = cc_ref(1, k, 2) + cc_ref(1, k, 4);
        tr2 = cc_ref(1, k, 1) + cc_ref(1, k, 3);
        ch_ref(1, 1, k) = tr1 + tr2;
        ch_ref(ido, 4, k) = tr2 - tr1;
        ch_ref(ido, 2, k) = cc_ref(1, k, 1) - cc_ref(1, k, 3);
        ch_ref(1, 3, k) = cc_ref(1, k, 4) - cc_ref(1, k, 2);
    }
    if (ido < 2) return;
    if (ido != 2) {
        idp2 = ido + 2;
        for (k = 1; k <= l1; ++k) {
            for (i = 3; i <= ido; i += 2) {
                ic = idp2 - i;
                cr2 = wa1[i - 2] * cc_ref(i - 1, k, 2) + wa1[i - 1] *
                                                         cc_ref(i, k, 2);
                ci2 = wa1[i - 2] * cc_ref(i, k, 2) - wa1[i - 1] * cc_ref(
                        i - 1, k, 2);
                cr3 = wa2[i - 2] * cc_ref(i - 1, k, 3) + wa2[i - 1] *
                                                         cc_ref(i, k, 3);
                ci3 = wa2[i - 2] * cc_ref(i, k, 3) - wa2[i - 1] * cc_ref(
                        i - 1, k, 3);
                cr4 = wa3[i - 2] * cc_ref(i - 1, k, 4) + wa3[i - 1] *
                                                         cc_ref(i, k, 4);
                ci4 = wa3[i - 2] * cc_ref(i, k, 4) - wa3[i - 1] * cc_ref(
                        i - 1, k, 4);
                tr1 = cr2 + cr4;
                tr4 = cr4 - cr2;
                ti1 = ci2 + ci4;
                ti4 = ci2 - ci4;
                ti2 = cc_ref(i, k, 1) + ci3;
                ti3 = cc_ref(i, k, 1) - ci3;
                tr2 = cc_ref(i - 1, k, 1) + cr3;
                tr3 = cc_ref(i - 1, k, 1) - cr3;
                ch_ref(i - 1, 1, k) = tr1 + tr2;
                ch_ref(ic - 1, 4, k) = tr2 - tr1;
                ch_ref(i, 1, k) = ti1 + ti2;
                ch_ref(ic, 4, k) = ti1 - ti2;
                ch_ref(i - 1, 3, k) = ti4 + tr3;
                ch_ref(ic - 1, 2, k) = tr3 - ti4;
                ch_ref(i, 3, k) = tr4 + ti3;
                ch_ref(ic, 2, k) = tr4 - ti3;
            }
        }
        if (ido % 2 == 1) {
            return;
        }
    }
    for (k = 1; k <= l1; ++k) {
        ti1 = -hsqt2 * (cc_ref(ido, k, 2) + cc_ref(ido, k, 4));
        tr1 = hsqt2 * (cc_ref(ido, k, 2) - cc_ref(ido, k, 4));
        ch_ref(ido, 1, k) = tr1 + cc_ref(ido, k, 1);
        ch_ref(ido, 3, k) = cc_ref(ido, k, 1) - tr1;
        ch_ref(1, 2, k) = ti1 - cc_ref(ido, k, 3);
        ch_ref(1, 4, k) = ti1 + cc_ref(ido, k, 3);
    }
}

#undef ch_ref
#undef cc_ref


static void radf5(int ido, int l1, const real *cc, real *ch,
                  const real *wa1, const real *wa2, const real *wa3, const real *wa4) {

    static const real tr11 = .309016994374947f;
    static const real ti11 = .951056516295154f;
    static const real tr12 = -.809016994374947f;
    static const real ti12 = .587785252292473f;

    int cc_offset, ch_offset;

    int i, k, ic;
    real ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3, dr4, dr5,
            cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
    int idp2;


#define cc_ref(a_1, a_2, a_3) cc[((a_3)*l1 + (a_2))*ido + a_1]
#define ch_ref(a_1, a_2, a_3) ch[((a_3)*5 + (a_2))*ido + a_1]

    ch_offset = 1 + ido * 6;
    ch -= ch_offset;
    cc_offset = 1 + ido * (1 + l1);
    cc -= cc_offset;
    --wa1;
    --wa2;
    --wa3;
    --wa4;

    for (k = 1; k <= l1; ++k) {
        cr2 = cc_ref(1, k, 5) + cc_ref(1, k, 2);
        ci5 = cc_ref(1, k, 5) - cc_ref(1, k, 2);
        cr3 = cc_ref(1, k, 4) + cc_ref(1, k, 3);
        ci4 = cc_ref(1, k, 4) - cc_ref(1, k, 3);
        ch_ref(1, 1, k) = cc_ref(1, k, 1) + cr2 + cr3;
        ch_ref(ido, 2, k) = cc_ref(1, k, 1) + tr11 * cr2 + tr12 * cr3;
        ch_ref(1, 3, k) = ti11 * ci5 + ti12 * ci4;
        ch_ref(ido, 4, k) = cc_ref(1, k, 1) + tr12 * cr2 + tr11 * cr3;
        ch_ref(1, 5, k) = ti12 * ci5 - ti11 * ci4;
    }
    if (ido == 1) {
        return;
    }
    idp2 = ido + 2;
    for (k = 1; k <= l1; ++k) {
        for (i = 3; i <= ido; i += 2) {
            ic = idp2 - i;
            dr2 = wa1[i - 2] * cc_ref(i - 1, k, 2) + wa1[i - 1] * cc_ref(i, k, 2);
            di2 = wa1[i - 2] * cc_ref(i, k, 2) - wa1[i - 1] * cc_ref(i - 1, k, 2);
            dr3 = wa2[i - 2] * cc_ref(i - 1, k, 3) + wa2[i - 1] * cc_ref(i, k, 3);
            di3 = wa2[i - 2] * cc_ref(i, k, 3) - wa2[i - 1] * cc_ref(i - 1, k, 3);
            dr4 = wa3[i - 2] * cc_ref(i - 1, k, 4) + wa3[i - 1] * cc_ref(i, k, 4);
            di4 = wa3[i - 2] * cc_ref(i, k, 4) - wa3[i - 1] * cc_ref(i - 1, k, 4);
            dr5 = wa4[i - 2] * cc_ref(i - 1, k, 5) + wa4[i - 1] * cc_ref(i, k, 5);
            di5 = wa4[i - 2] * cc_ref(i, k, 5) - wa4[i - 1] * cc_ref(i - 1, k, 5);
            cr2 = dr2 + dr5;
            ci5 = dr5 - dr2;
            cr5 = di2 - di5;
            ci2 = di2 + di5;
            cr3 = dr3 + dr4;
            ci4 = dr4 - dr3;
            cr4 = di3 - di4;
            ci3 = di3 + di4;
            ch_ref(i - 1, 1, k) = cc_ref(i - 1, k, 1) + cr2 + cr3;
            ch_ref(i, 1, k) = cc_ref(i, k, 1) + ci2 + ci3;
            tr2 = cc_ref(i - 1, k, 1) + tr11 * cr2 + tr12 * cr3;
            ti2 = cc_ref(i, k, 1) + tr11 * ci2 + tr12 * ci3;
            tr3 = cc_ref(i - 1, k, 1) + tr12 * cr2 + tr11 * cr3;
            ti3 = cc_ref(i, k, 1) + tr12 * ci2 + tr11 * ci3;
            tr5 = ti11 * cr5 + ti12 * cr4;
            ti5 = ti11 * ci5 + ti12 * ci4;
            tr4 = ti12 * cr5 - ti11 * cr4;
            ti4 = ti12 * ci5 - ti11 * ci4;
            ch_ref(i - 1, 3, k) = tr2 + tr5;
            ch_ref(ic - 1, 2, k) = tr2 - tr5;
            ch_ref(i, 3, k) = ti2 + ti5;
            ch_ref(ic, 2, k) = ti5 - ti2;
            ch_ref(i - 1, 5, k) = tr3 + tr4;
            ch_ref(ic - 1, 4, k) = tr3 - tr4;
            ch_ref(i, 5, k) = ti3 + ti4;
            ch_ref(ic, 4, k) = ti4 - ti3;
        }
    }
}

#undef ch_ref
#undef cc_ref


static void radfg(int ido, int ip, int l1, int idl1,
                  real *cc, real *c1, real *c2, real *ch, real *ch2, const real *wa) {
    int ch_offset, cc_offset,
            c1_offset, c2_offset, ch2_offset;

    int i, j, k, l, j2, ic, jc, lc, ik, is;
    real dc2, ai1, ai2, ar1, ar2, ds2;
    int nbd;
    real dcp, arg, dsp, ar1h, ar2h;
    int idp2, ipp2, idij, ipph;


#define c1_ref(a_1, a_2, a_3) c1[((a_3)*l1 + (a_2))*ido + a_1]
#define c2_ref(a_1, a_2) c2[(a_2)*idl1 + a_1]
#define cc_ref(a_1, a_2, a_3) cc[((a_3)*ip + (a_2))*ido + a_1]
#define ch_ref(a_1, a_2, a_3) ch[((a_3)*l1 + (a_2))*ido + a_1]
#define ch2_ref(a_1, a_2) ch2[(a_2)*idl1 + a_1]

    ch_offset = 1 + ido * (1 + l1);
    ch -= ch_offset;
    c1_offset = 1 + ido * (1 + l1);
    c1 -= c1_offset;
    cc_offset = 1 + ido * (1 + ip);
    cc -= cc_offset;
    ch2_offset = 1 + idl1;
    ch2 -= ch2_offset;
    c2_offset = 1 + idl1;
    c2 -= c2_offset;
    --wa;

    arg = (2 * M_PI) / (real) (ip);
    dcp = cos(arg);
    dsp = sin(arg);
    ipph = (ip + 1) / 2;
    ipp2 = ip + 2;
    idp2 = ido + 2;
    nbd = (ido - 1) / 2;
    if (ido == 1) {
        for (ik = 1; ik <= idl1; ++ik) {
            c2_ref(ik, 1) = ch2_ref(ik, 1);
        }
    } else {
        for (ik = 1; ik <= idl1; ++ik) {
            ch2_ref(ik, 1) = c2_ref(ik, 1);
        }
        for (j = 2; j <= ip; ++j) {
            for (k = 1; k <= l1; ++k) {
                ch_ref(1, k, j) = c1_ref(1, k, j);
            }
        }
        if (nbd <= l1) {
            is = -(ido);
            for (j = 2; j <= ip; ++j) {
                is += ido;
                idij = is;
                for (i = 3; i <= ido; i += 2) {
                    idij += 2;
                    for (k = 1; k <= l1; ++k) {
                        ch_ref(i - 1, k, j) = wa[idij - 1] * c1_ref(i - 1, k, j)
                                              + wa[idij] * c1_ref(i, k, j);
                        ch_ref(i, k, j) = wa[idij - 1] * c1_ref(i, k, j) - wa[idij] * c1_ref(i - 1, k, j);
                    }
                }
            }
        } else {
            is = -(ido);
            for (j = 2; j <= ip; ++j) {
                is += ido;
                for (k = 1; k <= l1; ++k) {
                    idij = is;
                    for (i = 3; i <= ido; i += 2) {
                        idij += 2;
                        ch_ref(i - 1, k, j) = wa[idij - 1] * c1_ref(i - 1, k, j) + wa[idij] * c1_ref(i, k, j);
                        ch_ref(i, k, j) = wa[idij - 1] * c1_ref(i, k, j) - wa[idij] * c1_ref(i - 1, k, j);
                    }
                }
            }
        }
        if (nbd >= l1) {
            for (j = 2; j <= ipph; ++j) {
                jc = ipp2 - j;
                for (k = 1; k <= l1; ++k) {
                    for (i = 3; i <= ido; i += 2) {
                        c1_ref(i - 1, k, j) = ch_ref(i - 1, k, j) + ch_ref(i -
                                                                           1, k, jc);
                        c1_ref(i - 1, k, jc) = ch_ref(i, k, j) - ch_ref(i, k,
                                                                        jc);
                        c1_ref(i, k, j) = ch_ref(i, k, j) + ch_ref(i, k, jc);
                        c1_ref(i, k, jc) = ch_ref(i - 1, k, jc) - ch_ref(i - 1,
                                                                         k, j);
                    }
                }
            }
        } else {
            for (j = 2; j <= ipph; ++j) {
                jc = ipp2 - j;
                for (i = 3; i <= ido; i += 2) {
                    for (k = 1; k <= l1; ++k) {
                        c1_ref(i - 1, k, j) = ch_ref(i - 1, k, j) + ch_ref(i -
                                                                           1, k, jc);
                        c1_ref(i - 1, k, jc) = ch_ref(i, k, j) - ch_ref(i, k,
                                                                        jc);
                        c1_ref(i, k, j) = ch_ref(i, k, j) + ch_ref(i, k, jc);
                        c1_ref(i, k, jc) = ch_ref(i - 1, k, jc) - ch_ref(i - 1,
                                                                         k, j);
                    }
                }
            }
        }
    }
    for (j = 2; j <= ipph; ++j) {
        jc = ipp2 - j;
        for (k = 1; k <= l1; ++k) {
            c1_ref(1, k, j) = ch_ref(1, k, j) + ch_ref(1, k, jc);
            c1_ref(1, k, jc) = ch_ref(1, k, jc) - ch_ref(1, k, j);
        }
    }

    ar1 = 1.f;
    ai1 = 0.f;
    for (l = 2; l <= ipph; ++l) {
        lc = ipp2 - l;
        ar1h = dcp * ar1 - dsp * ai1;
        ai1 = dcp * ai1 + dsp * ar1;
        ar1 = ar1h;
        for (ik = 1; ik <= idl1; ++ik) {
            ch2_ref(ik, l) = c2_ref(ik, 1) + ar1 * c2_ref(ik, 2);
            ch2_ref(ik, lc) = ai1 * c2_ref(ik, ip);
        }
        dc2 = ar1;
        ds2 = ai1;
        ar2 = ar1;
        ai2 = ai1;
        for (j = 3; j <= ipph; ++j) {
            jc = ipp2 - j;
            ar2h = dc2 * ar2 - ds2 * ai2;
            ai2 = dc2 * ai2 + ds2 * ar2;
            ar2 = ar2h;
            for (ik = 1; ik <= idl1; ++ik) {
                ch2_ref(ik, l) = ch2_ref(ik, l) + ar2 * c2_ref(ik, j);
                ch2_ref(ik, lc) = ch2_ref(ik, lc) + ai2 * c2_ref(ik, jc);
            }
        }
    }
    for (j = 2; j <= ipph; ++j) {
        for (ik = 1; ik <= idl1; ++ik) {
            ch2_ref(ik, 1) = ch2_ref(ik, 1) + c2_ref(ik, j);
        }
    }

    if (ido >= l1) {
        for (k = 1; k <= l1; ++k) {
            for (i = 1; i <= ido; ++i) {
                cc_ref(i, 1, k) = ch_ref(i, k, 1);
            }
        }
    } else {
        for (i = 1; i <= ido; ++i) {
            for (k = 1; k <= l1; ++k) {
                cc_ref(i, 1, k) = ch_ref(i, k, 1);
            }
        }
    }
    for (j = 2; j <= ipph; ++j) {
        jc = ipp2 - j;
        j2 = j + j;
        for (k = 1; k <= l1; ++k) {
            cc_ref(ido, j2 - 2, k) = ch_ref(1, k, j);
            cc_ref(1, j2 - 1, k) = ch_ref(1, k, jc);
        }
    }
    if (ido == 1) {
        return;
    }
    if (nbd >= l1) {
        for (j = 2; j <= ipph; ++j) {
            jc = ipp2 - j;
            j2 = j + j;
            for (k = 1; k <= l1; ++k) {
                for (i = 3; i <= ido; i += 2) {
                    ic = idp2 - i;
                    cc_ref(i - 1, j2 - 1, k) = ch_ref(i - 1, k, j) + ch_ref(
                            i - 1, k, jc);
                    cc_ref(ic - 1, j2 - 2, k) = ch_ref(i - 1, k, j) - ch_ref(
                            i - 1, k, jc);
                    cc_ref(i, j2 - 1, k) = ch_ref(i, k, j) + ch_ref(i, k, jc);
                    cc_ref(ic, j2 - 2, k) = ch_ref(i, k, jc) - ch_ref(i, k, j);
                }
            }
        }
    } else {
        for (j = 2; j <= ipph; ++j) {
            jc = ipp2 - j;
            j2 = j + j;
            for (i = 3; i <= ido; i += 2) {
                ic = idp2 - i;
                for (k = 1; k <= l1; ++k) {
                    cc_ref(i - 1, j2 - 1, k) = ch_ref(i - 1, k, j) + ch_ref(
                            i - 1, k, jc);
                    cc_ref(ic - 1, j2 - 2, k) = ch_ref(i - 1, k, j) - ch_ref(
                            i - 1, k, jc);
                    cc_ref(i, j2 - 1, k) = ch_ref(i, k, j) + ch_ref(i, k,
                                                                    jc);
                    cc_ref(ic, j2 - 2, k) = ch_ref(i, k, jc) - ch_ref(i, k, j);
                }
            }
        }
    }
}

#undef ch2_ref
#undef ch_ref
#undef cc_ref
#undef c2_ref
#undef c1_ref


static int decompose(int N, int *p) {
    int n = 0;

    while (N % 4 == 0) {
        N /= 4;
        p[n++] = 4;
    }

    while (N % 2 == 0) {
        N /= 2;
        p[n++] = 2;
    }
    for (int i = 3; i * i <= N; i += 2) {
        while (N % i == 0) {
            N /= i;
            p[n++] = i;
        }
    }
    if (N > 1) {
        p[n++] = N;
    }
    return n;
}

void rfftb(PACK_HEADER *header) {
    int i, k1, l1, l2, na, nf, ip, iw, ix2, ix3, ix4, ido, idl1;

    real *in = *(header->freq_domain);
    real *out = *(header->time_domain);
    nf = header->P_n;
    real *wa = header->W;
    int n = header->N;
    na = 0;
    l1 = 1;
    iw = 0;
    for (k1 = 1; k1 <= nf; ++k1) {
        ip = header->P[k1 - 1];
        l2 = ip * l1;
        ido = n / l2;
        idl1 = ido * l1;
        switch (ip) {
            case 4:
                ix2 = iw + ido;
                ix3 = ix2 + ido;
                radb4(ido, l1, na ? out : in, na ? in : out, &wa[iw], &wa[ix2], &wa[ix3]);
                na = 1 - na;
                break;
            case 2:
                radb2(ido, l1, na ? out : in, na ? in : out, &wa[iw]);
                na = 1 - na;
                break;
            case 3:
                ix2 = iw + ido;
                radb3(ido, l1, na ? out : in, na ? in : out, &wa[iw], &wa[ix2]);
                na = 1 - na;
                break;
            case 5:
                ix2 = iw + ido;
                ix3 = ix2 + ido;
                ix4 = ix3 + ido;
                radb5(ido, l1, na ? out : in, na ? in : out, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
                na = 1 - na;
                break;
            default:
                if (na == 0) {
                    radbg(ido, ip, l1, idl1, in, in, in, out, out, &wa[iw]);
                } else {
                    radbg(ido, ip, l1, idl1, out, out, out, in, in, &wa[iw]);
                }
                if (ido == 1) {
                    na = 1 - na;
                }
                break;
        }
        l1 = l2;
        iw += (ip - 1) * ido;
    }
    if (na == 0) {

        *(header->time_domain) = in;
        *(header->freq_domain) = out;

    }
}

void rfftf(PACK_HEADER *header) {
    int i, k1, l1, l2, na, kh, nf, ip, iw, ix2, ix3, ix4, ido, idl1;

    real *in = *(header->time_domain);
    real *out = *(header->freq_domain);
    real *W = header->W;
    int N = header->N;
    nf = header->P_n;
    na = 1;
    l2 = N;
    iw = N - 1;
    for (k1 = 1; k1 <= nf; ++k1) {
        kh = nf - k1;
        ip = header->P[kh];
        l1 = l2 / ip;
        ido = N / l2;
        idl1 = ido * l1;
        iw -= (ip - 1) * ido;
        na = 1 - na;
        switch (ip) {
            case 4:
                ix2 = iw + ido;
                ix3 = ix2 + ido;
                radf4(ido, l1, na ? out : in, na ? in : out, &W[iw], &W[ix2], &W[ix3]);
                break;
            case 2:
                radf2(ido, l1, na ? out : in, na ? in : out, &W[iw]);
                break;
            case 3:
                ix2 = iw + ido;
                radf3(ido, l1, na ? out : in, na ? in : out, &W[iw], &W[ix2]);
                break;
            case 5:
                ix2 = iw + ido;
                ix3 = ix2 + ido;
                ix4 = ix3 + ido;
                radf5(ido, l1, na ? out : in, na ? in : out, &W[iw], &W[ix2], &W[ix3], &W[ix4]);
                break;
            default:
                if (ido == 1) {
                    na = 1 - na;
                }
                if (na == 0) {
                    radfg(ido, ip, l1, idl1, in, in, in, out, out, &W[iw]);
                    na = 1;
                } else {
                    radfg(ido, ip, l1, idl1, out, out, out, in, in, &W[iw]);
                    na = 0;
                }
                break;
        }
        l2 = l1;
    }
    if (na == 1) {
        *(header->freq_domain) = in;
        *(header->time_domain) = out;
    }

}

PACK_HEADER *pack_new(real **time_domain, real **freq_domain, int N) {
    PACK_HEADER *p = (PACK_HEADER *) malloc(sizeof(PACK_HEADER));
    p->N = N;
    p->P_n = decompose(N, p->P);
    p->W = (real *) malloc(sizeof(real) * N);
    p->time_domain = time_domain;
    p->freq_domain = freq_domain;

    int i, l1, l2;
    real fi;
    int ld, ii, nf, ip, is;
    real arg;
    int ido, ipm;
    int nfm1;
    real argh;
    real arg_ld;


    argh = (2 * M_PI) / (real) (p->N);
    is = 0;
    nfm1 = p->P_n - 1;
    l1 = 1;
    if (nfm1 == 0) {
        return p;
    }
    for (int k = 0; k < nfm1; ++k) {
        ip = p->P[k];
        ld = 0;
        l2 = l1 * ip;
        ido = p->N / l2;
        ipm = ip - 1;
        for (int j = 0; j < ipm; ++j) {
            ld += l1;
            i = is;
            arg_ld = (real) ld * argh;
            fi = 0.f;
            for (ii = 3; ii <= ido; ii += 2) {
                i += 2;
                fi += 1.f;
                arg = fi * arg_ld;
                p->W[i - 2] = real_cos(arg);
                p->W[i - 1] = real_sin(arg);

            }
            is += ido;
        }
        l1 = l2;
    }
    return p;
}

void pack_fft(PACK_HEADER *header, fft_type type) {
    switch (type) {
        case Forward:
            rfftf(header);
            break;
        case Inverse:
            rfftb(header);
            break;
    }
}

void fft_pack_real_test() {
    int N = 15625;


    real *time_domain = (real *) malloc(N * sizeof(real));
    real *freq_domain = (real *) malloc(N * sizeof(real));


    for (int i = 0; i < N / 2; i++) {
        time_domain[i] = i + 1 % 51;
        freq_domain[i] = 0;
    }
    for (int i = N / 2; i < N; ++i) {
        time_domain[i] = 0;
        freq_domain[i] = 0;
    }


    PACK_HEADER *p = pack_new(&time_domain, &freq_domain, N);

//    // start timer
//    clock_t start = clock();



    pack_fft(p, Forward);

//
//    freq_domain[0] = freq_domain[0] * freq_domain[0];
//    freq_domain[N - 1] = freq_domain[N - 1] * freq_domain[N - 1];
//    real c, d;
//    for (int i = 1; i < (N) / 2; ++i) {
//        c = freq_domain[2 * i - 1] * freq_domain[2 * i - 1] - freq_domain[2 * i] * freq_domain[2 * i];
//        d = freq_domain[2 * i - 1] * freq_domain[2 * i] + freq_domain[2 * i - 1] * freq_domain[2 * i];
//        freq_domain[2 * i - 1] = c;
//        freq_domain[2 * i] = d;
//    }



    pack_fft(p, Inverse);

    // stop timer
//    clock_t end = clock();
//


//    for (int i = 0; i < N % 51; ++i) {
//        printf("%f \n", time_domain[i] / N);
//    }
}

void pack_free(PACK_HEADER *header) {
    free(header->W);
    free(header);
}
