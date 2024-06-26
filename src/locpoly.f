c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine locpol.f cccccccccc

c For computing an binned approximation to a
c local bandwidth local polynomial kernel regression estimator
c of an arbitrary derivative of a regression function.
c LINPACK is used for matrix inversion.

c Last changed: 10/02/95

      subroutine locpol(xcnts,ycnts,idrv,delta,hdisc,Lvec,indic,
     +                  midpts,M,iQ,fkap,ipp,ippp,ss,tt,Smat,Tvec,
     +                  ipvt,cvest)
      integer i,j,k,ii,Lvec(*),M,iQ,mid,indic(*),midpts(*),ipvt(*),
     +        info,idrv,ipp,ippp,indss
      double precision xcnts(*),ycnts(*),fkap(*),hdisc(*),
     +                 cvest(*),delta,ss(M,ippp),tt(M,ipp),
     +                 Smat(ipp,ipp),Tvec(ipp),fac

c Obtain kernel weights

      mid = Lvec(1) + 1
      do i=1,(iQ-1)
         midpts(i) = mid
         fkap(mid) = 1.0d0
         do j=1,Lvec(i)
            fkap(mid+j) = exp(-(delta*j/hdisc(i))**2/2)
            fkap(mid-j) = fkap(mid+j)
         end do
         mid = mid + Lvec(i) + Lvec(i+1) + 1
      end do
      midpts(iQ) = mid
      fkap(mid) = 1.0d0
      do j=1,Lvec(iQ)
         fkap(mid+j) = exp(-(delta*j/hdisc(iQ))**2/2)
         fkap(mid-j) = fkap(mid+j)
      end do

c Combine kernel weights and grid counts

      do k = 1,M
         if (xcnts(k).ne.0) then
            do i = 1,iQ
               do j = max(1,k-Lvec(i)),min(M,k+Lvec(i))
                  if (indic(j).eq.i) then
                     fac = 1.0d0
                     ss(j,1) = ss(j,1) + xcnts(k)*fkap(k-j+midpts(i))
                     tt(j,1) = tt(j,1) + ycnts(k)*fkap(k-j+midpts(i))
                     do ii = 2,ippp
                        fac = fac*delta*(k-j)
                        ss(j,ii) = ss(j,ii)
     +                  + xcnts(k)*fkap(k-j+midpts(i))*fac
                        if (ii.le.ipp) then
                           tt(j,ii) = tt(j,ii)
     +                     + ycnts(k)*fkap(k-j+midpts(i))*fac
                        endif
                     end do
                  endif
               end do
            end do
         endif
         end do

      do k = 1,M
         do i = 1,ipp
            do j = 1,ipp
               indss = i + j - 1
               Smat(i,j) = ss(k,indss)
            end do
            Tvec(i) = tt(k,i)
         end do

        call dgefa(Smat,ipp,ipp,ipvt,info)
        call dgesl(Smat,ipp,ipp,ipvt,Tvec,0)

        cvest(k) = Tvec(idrv+1)

      end do

      return
      end

cccccccccc End of locpol.f cccccccccc
