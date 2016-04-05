#include "ml_ndt_scanmatching/other/layer.h"

//namespace ml_ndt{

  // perform line search to find the best descent rate (More&Thuente)
  double Layer::lineSearchMT(const transform_t & trans,
                              const pose_t & gradient,
                              const pose_t & inc,
                              const points_t & cloud_in) const {
    // default params
    double stp = 4.0; // default step
    double recoverystep = 0.05;
    double dginit = 0.0;
    double ftol = 0.0001; // epsilon 1
    double gtol = 0.9999; // epsilon 2
    double stpmax = 10.0;
    double stpmin = 0.0001;
    int maxfev = 10; // max function evaluations
    double xtol = 0.01; // window of uncertainty around the optimal step

    double direction = 1.0;
    double score_init = 0.0;

    Eigen::Transform<double,2, Eigen::TransformTraits::Affine> ps, ps2;
    pose_t pincr, score_gradient_here,increment;
    /////
    increment = inc;
    int info = 0;  // return code
    int infoc = 1; // return code for subroutine cstep

    // Compute the initial gradient in the search direction and check
    // that s is a descent direction.

    // we want to maximize s, so we should minimize -s
    score_init = scoreLayer(trans, cloud_in);

    // gradient directions are opposite for the negated function
    // gradient = -gradient;

    dginit = increment.dot(gradient);
    if (dginit >= 0.0) {
      //    cout << "MoreThuente::cvsrch - wrong direction (dginit = " << dginit
      //    << ")" << endl;
      // return recoverystep; //TODO TSV -1; //

      increment = -increment;
      dginit = -dginit;
      direction = -1;

      if (dginit >= 0.0) {
        //    cout << "MoreThuente::cvsrch - Non-descent direction (dginit = " <<
        //    dginit << ")" << endl;
        // stp = recoverystep;
        // newgrp.computeX(oldgrp, dir, stp);
        return recoverystep;
      }
    } else {
      //     cout<<"correct direction (dginit = " << dginit << ")" << endl;
    }

    // Initialize local variables.

    bool brackt = false;            // has the soln been bracketed?
    bool stage1 = true;             // are we in stage 1?
    int nfev = 0;                   // number of function evaluations
    double dgtest = ftol * dginit;  // f for curvature condition
    double width = stpmax - stpmin; // interval width
    double width1 = 2 * width;      // ???

    // cout<<"dgtest "<<dgtest<<endl;
    // initial function value
    double finit = 0.0;
    finit = score_init;

    // The variables stx, fx, dgx contain the values of the step,
    // function, and directional derivative at the best step.  The
    // variables sty, fy, dgy contain the value of the step, function,
    // and derivative at the other endpoint of the interval of
    // uncertainty.  The variables stp, f, dg contain the values of the
    // step, function, and derivative at the current step.

    double stx = 0.0;
    double fx = finit;
    double dgx = dginit;
    double sty = 0.0;
    double fy = finit;
    double dgy = dginit;

    // Get the linear solve tolerance for adjustable forcing term
    double eta_original = -1.0;
    double eta = 0.0;
    eta = eta_original;

    // Start of iteration.

    double stmin, stmax;
    double fm, fxm, fym, dgm, dgxm, dgym;

    while (1) {
      // Set the minimum and maximum steps to correspond to the present
      // interval of uncertainty.
      if (brackt) {
        stmin = MoreThuente::min(stx, sty);
        stmax = MoreThuente::max(stx, sty);
      } else {
        stmin = stx;
        stmax = stp + 4 * (stp - stx);
      }

      // Force the step to be within the bounds stpmax and stpmin.
      stp = MoreThuente::max(stp, stpmin);
      stp = MoreThuente::min(stp, stpmax);

      // If an unusual termination is to occur then let stp be the
      // lowest point obtained so far.

      if ((brackt && ((stp <= stmin) || (stp >= stmax))) ||
          (nfev >= maxfev - 1) || (infoc == 0) ||
          (brackt && (stmax - stmin <= xtol * stmax))) {
        stp = stx;
      }

      // Evaluate the function and gradient at stp
      // and compute the directional derivative.
      ///////////////////////////////////////////////////////////////////////////

      pincr = stp * increment;

      ps = ps.fromPositionOrientationScale(pincr.head<2>(),
                                            Eigen::Rotation2Dd(pincr(2)),
                                            Eigen::Vector2d::Ones());
      double f = 0.0;
      double total_score=0.0;
      Eigen::Vector3d g = Eigen::Vector3d::Zero();
      Eigen::Matrix<double, 2, 3> jacobian;
  
      Eigen::Rotation2D<double> rot(0);
      rot = rot.fromRotationMatrix(ps.rotation());
      double si = std::sin(rot.angle());
      double co = std::cos(rot.angle());
      score_gradient_here.setZero();
      ps2.setIdentity();

      for (const auto & point : cloud_in) {
        Eigen::Vector2d trans_point = ps*point;
        double x = trans_point(0);
        double y = trans_point(1);
        Field field;
        if(!getPointField(trans_point,field))
          continue;
        Eigen::Matrix2d inv_covar =field.getInvCovar();
        //DEBUG(inv_variace);
        point_t difference = trans_point - field.getMean();
        double point_score = scorePoint(field,trans_point);
        //DEBUG(difference.dot(field.calcInvertedVariance() * difference)* -0.5F);
        //DEBUG(point_score);
        total_score += point_score;
        jacobian << 1, 0, -x * si - y * co, 
                    0, 1,  x * co - y * si;
        g += pointGradient(difference,inv_covar,point_score,jacobian);
      }
      f = total_score;
      score_gradient_here = g;

      // VALGRIND_CHECK_VALUE_IS_DEFINED(score_gradient_here);
      // VALGRIND_CHECK_VALUE_IS_DEFINED(increment);
      double dg = 0.0;
      dg = increment.dot(score_gradient_here);

      // VALGRIND_CHECK_VALUE_IS_DEFINED(dg);
      // cout<<"dg = "<<dg<<endl;
      nfev++;

      ///////////////////////////////////////////////////////////////////////////

      // cout<<"consider step "<<stp<<endl;
      // Armijo-Goldstein sufficient decrease
      double ftest1 = finit + stp * dgtest;
      // cout<<"ftest1 is "<<ftest1<<endl;

      // Test for convergence.

      if ((brackt && ((stp <= stmin) || (stp >= stmax))) || (infoc == 0))
        info = 6; // Rounding errors

      if ((stp == stpmax) && (f <= ftest1) && (dg <= dgtest))
        info = 5; // stp=stpmax

      if ((stp == stpmin) && ((f > ftest1) || (dg >= dgtest)))
        info = 4; // stp=stpmin

      if (nfev >= maxfev)
        info = 3; // max'd out on fevals

      if (brackt && (stmax - stmin <= xtol * stmax))
        info = 2; // bracketed soln

      // RPP sufficient decrease test can be different
      bool sufficientDecreaseTest = false;
      sufficientDecreaseTest = (f <= ftest1); // Armijo-Golstein

      // cout<<"ftest2 "<<gtol*(-dginit)<<endl;
      // cout<<"sufficientDecrease? "<<sufficientDecreaseTest<<endl;
      // cout<<"curvature ok? "<<(fabs(dg) <= gtol*(-dginit))<<endl;
      if ((sufficientDecreaseTest) && (fabs(dg) <= gtol * (-dginit)))
        info = 1; // Success!!!!

      if (info != 0) // Line search is done
      {
        if (info != 1) // Line search failed
        {
          // RPP add
          // counter.incrementNumFailedLineSearches();

          // if (recoveryStepType == Constant)
          stp = recoverystep;

          // newgrp.computeX(oldgrp, dir, stp);

          // message = "(USING RECOVERY STEP!)";

        } else // Line search succeeded
        {
          // message = "(STEP ACCEPTED!)";
        }

        // print.printStep(nfev, stp, finit, f, message);

        // Returning the line search flag
        // cout<<"LineSearch::"<<message<<" info "<<info<<endl;
        return stp;

      } // info != 0

      // RPP add
      // counter.incrementNumIterations();

      // In the first stage we seek a step for which the modified
      // function has a nonpositive value and nonnegative derivative.

      if (stage1 && (f <= ftest1) &&
          (dg >= MoreThuente::min(ftol, gtol) * dginit)) {
        stage1 = false;
      }

      // A modified function is used to predict the step only if we have
      // not obtained a step for which the modified function has a
      // nonpositive function value and nonnegative derivative, and if a
      // lower function value has been obtained but the decrease is not
      // sufficient.

      if (stage1 && (f <= fx) && (f > ftest1)) {

        // Define the modified function and derivative values.

        fm = f - stp * dgtest;
        fxm = fx - stx * dgtest;
        fym = fy - sty * dgtest;
        dgm = dg - dgtest;
        dgxm = dgx - dgtest;
        dgym = dgy - dgtest;

        // Call cstep to update the interval of uncertainty
        // and to compute the new step.

        // VALGRIND_CHECK_VALUE_IS_DEFINED(dgm);
        infoc = MoreThuente::cstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm,
                                   brackt, stmin, stmax);

        // Reset the function and gradient values for f.

        fx = fxm + stx * dgtest;
        fy = fym + sty * dgtest;
        dgx = dgxm + dgtest;
        dgy = dgym + dgtest;

      }

      else {

        // Call cstep to update the interval of uncertainty
        // and to compute the new step.

        // VALGRIND_CHECK_VALUE_IS_DEFINED(dg);
        infoc = MoreThuente::cstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt,
                                   stmin, stmax);
      }

      // Force a sufficient decrease in the size of the
      // interval of uncertainty.

      if (brackt) {
        if (fabs(sty - stx) >= 0.66 * width1)
          stp = stx + 0.5 * (sty - stx);
        width1 = width;
        width = fabs(sty - stx);
      }

    } // while-loop
  }

  int Layer::MoreThuente::cstep(double &stx, double &fx, double &dx,
                                        double &sty, double &fy, double &dy,
                                        double &stp, double &fp, double &dp,
                                        bool &brackt, double stmin,
                                        double stmax) {
    int info = 0;

    // Check the input parameters for errors.

    if ((brackt && ((stp <= MoreThuente::min(stx, sty)) ||
                    (stp >= MoreThuente::max(stx, sty)))) ||
        (dx * (stp - stx) >= 0.0) || (stmax < stmin))
      return info;

    // Determine if the derivatives have opposite sign.

    double sgnd = dp * (dx / fabs(dx));

    // First case. A higher function value.  The minimum is
    // bracketed. If the cubic step is closer to stx than the quadratic
    // step, the cubic step is taken, else the average of the cubic and
    // quadratic steps is taken.

    bool bound;
    double theta;
    double s;
    double gamma;
    double p, q, r;
    double stpc, stpq, stpf;

    if (fp > fx) {
      info = 1;
      bound = 1;
      theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      // VALGRIND_CHECK_VALUE_IS_DEFINED(theta);
      // VALGRIND_CHECK_VALUE_IS_DEFINED(dx);
      // VALGRIND_CHECK_VALUE_IS_DEFINED(dp);
      s = MoreThuente::absmax(theta, dx, dp);
      gamma = s * sqrt(((theta / s) * (theta / s)) - (dx / s) * (dp / s));
      if (stp < stx)
        gamma = -gamma;

      p = (gamma - dx) + theta;
      q = ((gamma - dx) + gamma) + dp;
      r = p / q;
      stpc = stx + r * (stp - stx);
      stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2) * (stp - stx);
      if (fabs(stpc - stx) < fabs(stpq - stx))
        stpf = stpc;
      else
        stpf = stpc + (stpq - stpc) / 2;

      brackt = true;
    }

    // Second case. A lower function value and derivatives of opposite
    // sign. The minimum is bracketed. If the cubic step is closer to
    // stx than the quadratic (secant) step, the cubic step is taken,
    // else the quadratic step is taken.

    else if (sgnd < 0.0) {
      info = 2;
      bound = false;
      theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      s = MoreThuente::absmax(theta, dx, dp);
      gamma = s * sqrt(((theta / s) * (theta / s)) - (dx / s) * (dp / s));
      if (stp > stx)
        gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dx;
      r = p / q;
      stpc = stp + r * (stx - stp);
      stpq = stp + (dp / (dp - dx)) * (stx - stp);
      if (fabs(stpc - stp) > fabs(stpq - stp))
        stpf = stpc;
      else
        stpf = stpq;
      brackt = true;
    }

    // Third case. A lower function value, derivatives of the same sign,
    // and the magnitude of the derivative decreases.  The cubic step is
    // only used if the cubic tends to infinity in the direction of the
    // step or if the minimum of the cubic is beyond stp. Otherwise the
    // cubic step is defined to be either stmin or stmax. The
    // quadratic (secant) step is also computed and if the minimum is
    // bracketed then the the step closest to stx is taken, else the
    // step farthest away is taken.

    else if (fabs(dp) < fabs(dx)) {
      info = 3;
      bound = true;
      theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      s = MoreThuente::absmax(theta, dx, dp);

      // The case gamma = 0 only arises if the cubic does not tend
      // to infinity in the direction of the step.

      gamma = s * sqrt(max(0, (theta / s) * (theta / s) - (dx / s) * (dp / s)));
      if (stp > stx)
        gamma = -gamma;

      p = (gamma - dp) + theta;
      q = (gamma + (dx - dp)) + gamma;
      r = p / q;
      if ((r < 0.0) && (gamma != 0.0))
        stpc = stp + r * (stx - stp);
      else if (stp > stx)
        stpc = stmax;
      else
        stpc = stmin;

      stpq = stp + (dp / (dp - dx)) * (stx - stp);
      if (brackt) {
        if (fabs(stp - stpc) < fabs(stp - stpq))
          stpf = stpc;
        else
          stpf = stpq;
      } else {
        if (fabs(stp - stpc) > fabs(stp - stpq))
          stpf = stpc;
        else
          stpf = stpq;
      }
    }

    // Fourth case. A lower function value, derivatives of the same
    // sign, and the magnitude of the derivative does not decrease. If
    // the minimum is not bracketed, the step is either stmin or
    // stmax, else the cubic step is taken.

    else {
      info = 4;
      bound = false;
      if (brackt) {
        theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
        s = MoreThuente::absmax(theta, dy, dp);
        gamma = s * sqrt(((theta / s) * (theta / s)) - (dy / s) * (dp / s));
        if (stp > sty)
          gamma = -gamma;
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + dy;
        r = p / q;
        stpc = stp + r * (sty - stp);
        stpf = stpc;
      } else if (stp > stx)
        stpf = stmax;
      else
        stpf = stmin;
    }

    // Update the interval of uncertainty. This update does not depend
    // on the new step or the case analysis above.

    if (fp > fx) {
      sty = stp;
      fy = fp;
      dy = dp;
    } else {
      if (sgnd < 0.0) {
        sty = stx;
        fy = fx;
        dy = dx;
      }
      stx = stp;
      fx = fp;
      dx = dp;
    }

    // Compute the new step and safeguard it.

    stpf = MoreThuente::min(stmax, stpf);
    stpf = MoreThuente::max(stmin, stpf);
    stp = stpf;
    if (brackt && bound) {
      if (sty > stx)
        stp = min(stx + 0.66 * (sty - stx), stp);
      else
        stp = max(stx + 0.66 * (sty - stx), stp);
    }

    return info;
  }

  double Layer::MoreThuente::min(double a, double b) {
    return (a < b ? a : b);
  }

  double Layer::MoreThuente::max(double a, double b) {
    return (a > b ? a : b);
  }

  double Layer::MoreThuente::absmax(double a, double b, double c) {
    a = fabs(a);
    b = fabs(b);
    c = fabs(c);

    if (a > b)
      return (a > c) ? a : c;
    else
      return (b > c) ? b : c;
  }

//}
