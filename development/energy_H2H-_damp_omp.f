c     vectorial polarization and partial charges (+,-2,+) on atoms and center
c
c     H- assumed to be first
      
      ! JM Changed: Nmax to Natoms
      function func(Natoms,P)
      implicit double precision(a-h,o-z)
      integer Natoms
      dimension P(3*Natoms)
      dimension x(Natoms),y(Natoms),z(Natoms)
      dimension polar(Natoms),q(Natoms)
      
      dimension Ex(Natoms),Ey(Natoms),Ez(Natoms)

      parameter(rdmp=9.D0,rcut=3.2D0)

      parameter(facrep=1000.D0)

c     parameters for all-atom Rydberg potential
      parameter(bohr=0.529177D0)
      parameter(q_H=0.4829d0,rH2=0.74796d0/bohr,VH2=0.177095d0)
      parameter(aH2=2.159129d0,deltaH2=1.57253d-2)
      parameter(sigma2=(2.958D0/bohr)**2,epsilon=36.7D0*3.165D-6)

!     JM Changed: parameters for H2-H^- Buckingham-type potential initialized here
      parameter(A=1.0305769968240777,b=0.86775454773152316)
      parameter(C6=184.32635427872907,C8=0.11904965708491388)
      parameter(polar_H2=0.787d0/bohr**3,polar_Hminus=35.5d0)

      dimension xH2(Natoms/2,3),yH2(Natoms/2,3),zH2(Natoms/2,3),qH2(3)

      ! JM changed: Initialize the coordinates and the polarization terms
      do i=1,Natoms
         x(i)=P(3*i-2)
         y(i)=P(3*i-1)
         z(i)=P(3*i)
         polar(i)=polar_H2
      enddo
      polar(1)=polar_Hminus
      
c     polarization
      
      erep=0.d0
      edisp=0.d0
      epol=0.D0
      eh2=0.D0
      ecoul=0.D0
      eintra=0.d0

      xHm=x(1)
      yHm=y(1)
      zHm=z(1)
      q(1)=-1.D0

c     change H2 coordinates

      nH2=(Natoms-1)/2

      do i=1,nH2
         xH2(i,1)=x(2*i)
         yH2(i,1)=y(2*i)
         zH2(i,1)=z(2*i)

         xH2(i,2)=x(2*i+1)
         yH2(i,2)=y(2*i+1)
         zH2(i,2)=z(2*i+1)

         xH2(i,3)=(xH2(i,1)+xH2(i,2))/2.D0
         yH2(i,3)=(yH2(i,1)+yH2(i,2))/2.D0
         zH2(i,3)=(zH2(i,1)+zH2(i,2))/2.D0
      enddo
      qH2(1)=q_H
      qH2(2)=q_H
      qH2(3)=-2.D0*q_H

c     intramolecular contribution

      do i=1,nH2
         xij=xH2(i,1)-xH2(i,2)
         yij=yH2(i,1)-yH2(i,2)
         zij=zH2(i,1)-zH2(i,2)

         r12_2=xij**2+yij**2+zij**2
         r12=dsqrt(r12_2)
         dum=aH2*(r12/rH2-1.D0)
         phi=(1.D0+dum+deltaH2*dum**3)
         v_intra=VH2-VH2*phi*exp(-dum)

         eintra=eintra+v_intra
      enddo

c     electric field

      do i=1,nH2+1
         Ex(i)=0.D0
         Ey(i)=0.D0
         Ez(i)=0.D0
      enddo

c     pair energies (+ gradient)

      ih=1
      
      do i=1,nH2
         do k=1,3
            xih=xhm-xH2(i,k)
            yih=yhm-yH2(i,k)
            zih=zhm-zH2(i,k)
            xih2=xih**2
            yih2=yih**2
            zih2=zih**2
            rih2=xih2+yih2+zih2
            rih=sqrt(rih2)
            rih3=rih2*rih

            ddamp=0.D0
            if (rih.lt.rcut) then
               damp=exp(-(1.D0-rcut/rih)**2)
               ddamp=2.D0*rcut*(1.D0-rcut/rih)
            else
               damp=1.D0
            endif
            
            ecoul=ecoul+damp*q(1)*qH2(k)/rih

            ecoul=ecoul+facrep/rih**12

c     field on ion

            Ex(1)=Ex(1)+damp*qH2(k)*xih/rih3
            Ey(1)=Ey(1)+damp*qH2(k)*yih/rih3
            Ez(1)=Ez(1)+damp*qH2(k)*zih/rih3
            
c     bare electric fields on H2 centers (site k=3)

            if (k.eq.3) then
               Ex(i+1)=Ex(i+1)-q(1)*xih/rih3
               Ey(i+1)=Ey(i+1)-q(1)*yih/rih3
               Ez(i+1)=Ez(i+1)-q(1)*zih/rih3
               
c     repulsion-dispersion (also if site k=3)

               erep=erep+A*exp(-b*rih)

               hcdmp=1.d0
               dhcdmp=0.d0
               if (rih.lt.rdmp) then
                  hcdmp=exp(-(rdmp/rih-1.d0)**2)
               endif
                  
               hcdisp=C6/rih2**3+C8/rih2**4
                  
               edisp=edisp-hcdisp*hcdmp
                  
            endif

         enddo
      enddo

c     H2-H2 interaction

      do i=1,nH2
         do j=i+1,nH2
            do k=1,3
               do l=1,3
                  xik_jl=xH2(i,k)-xH2(j,l)
                  yik_jl=yH2(i,k)-yH2(j,l)
                  zik_jl=zH2(i,k)-zH2(j,l)
                  rik_jl_2=xik_jl**2+yik_jl**2+zik_jl**2
                  rik_jl=sqrt(rik_jl_2)
                  dumq=qH2(k)*qH2(l)/rik_jl

                  ecoul=ecoul+dumq

                  dumq=-dumq/rik_jl_2

                  if (k.eq.3.and.l.eq.3) then ! LJ between c.o.m.

                     rik_jl_6=(sigma2/rik_jl_2)**3
                     eH2=eH2+4.D0*epsilon*rik_jl_6*(rik_jl_6-1.D0)

                  endif

               enddo
            enddo
         enddo
      enddo

c     polarization

      do i=1,nH2+1
         epol=epol-polar(i)*(Ex(i)**2+Ey(i)**2+Ez(i)**2)/2.D0
      enddo

      epot=eintra+ecoul+erep+edisp+eH2+epol

      func=epot
      
      return
      end

c----------------------------------------------------------------------------

      ! JM Changed: Nmax to Natoms
      subroutine dfunc(Natoms,P,dP)
      implicit double precision(a-h,o-z)
      integer Natoms
      dimension P(3*Natoms),dP(3*Natoms)
      dimension x(Natoms),y(Natoms),z(Natoms)
      dimension gx(Natoms),gy(Natoms),gz(Natoms)
      dimension gxpol(Natoms),gypol(Natoms),gzpol(Natoms)
      dimension polar(Natoms),q(Natoms)

      dimension Ex(Natoms),Ey(Natoms),Ez(Natoms)
      dimension dexx(3*Natoms,3*Natoms),dexy(3*Natoms,3*Natoms)
      dimension dexz(3*Natoms,3*Natoms)
      dimension deyy(3*Natoms,3*Natoms),deyz(3*Natoms,3*Natoms)
      dimension dezz(3*Natoms,3*Natoms)

      parameter(rdmp=9.D0,rcut=3.2D0)

      parameter(facrep=1000.D0)

c     parameters for all-atom Rydberg potential
      parameter(bohr=0.529177D0)
      parameter(q_H=0.4829d0,rH2=0.74796d0/bohr,VH2=0.177095d0)
      parameter(aH2=2.159129d0,deltaH2=1.57253d-2)
      parameter(sigma2=(2.958D0/bohr)**2,epsilon=36.7D0*3.165D-6)

!     JM changed: parameters for H2-H^- Buckingham-type potential initialized here
      parameter(A=1.0305769968240777,b=0.86775454773152316)
      parameter(C6=184.32635427872907,C8=0.11904965708491388)
      parameter(polar_H2=0.787d0/bohr**3,polar_Hminus=35.5d0)

      dimension xH2(Natoms/2,3),yH2(Natoms/2,3),zH2(Natoms/2,3),qH2(3)

      ! JM changed: Initialize the coordinates and the polarization terms
      do i=1,Natoms
         x(i)=P(3*i-2)
         y(i)=P(3*i-1)
         z(i)=P(3*i)
         polar(i)=polar_H2
      enddo
      polar(1)=polar_Hminus
      
c     polarization
      
      erep=0.d0
      edisp=0.d0
      epol=0.D0
      eh2=0.D0
      ecoul=0.D0
      eintra=0.d0

      xHm=x(1)
      yHm=y(1)
      zHm=z(1)
      q(1)=-1.D0

c     change H2 coordinates

      nH2=(Natoms-1)/2

      do i=1,nH2
         xH2(i,1)=x(2*i)
         yH2(i,1)=y(2*i)
         zH2(i,1)=z(2*i)

         xH2(i,2)=x(2*i+1)
         yH2(i,2)=y(2*i+1)
         zH2(i,2)=z(2*i+1)

         xH2(i,3)=(xH2(i,1)+xH2(i,2))/2.D0
         yH2(i,3)=(yH2(i,1)+yH2(i,2))/2.D0
         zH2(i,3)=(zH2(i,1)+zH2(i,2))/2.D0
      enddo
      qH2(1)=q_H
      qH2(2)=q_H
      qH2(3)=-2.D0*q_H

         do i=1,Natoms
            gx(i)=0.D0
            gy(i)=0.D0
            gz(i)=0.D0

            gxpol(i)=0.D0
            gypol(i)=0.D0
            gzpol(i)=0.D0
         enddo

c     intramolecular contribution

      do i=1,nH2
         xij=xH2(i,1)-xH2(i,2)
         yij=yH2(i,1)-yH2(i,2)
         zij=zH2(i,1)-zH2(i,2)

         r12_2=xij**2+yij**2+zij**2
         r12=dsqrt(r12_2)
         dum=aH2*(r12/rH2-1.D0)
         phi=(1.D0+dum+deltaH2*dum**3)
         v_intra=VH2-VH2*phi*exp(-dum)

         eintra=eintra+v_intra

            dphi=dum-3.D0*deltaH2*dum**2+deltaH2*dum**3
            dphi=VH2*dphi*(aH2/rH2)*exp(-dum)

            gx(2*i  )=gx(2*i  )+xij*dphi/r12
            gy(2*i  )=gy(2*i  )+yij*dphi/r12
            gz(2*i  )=gz(2*i  )+zij*dphi/r12

            gx(2*i+1)=gx(2*i+1)-xij*dphi/r12
            gy(2*i+1)=gy(2*i+1)-yij*dphi/r12
            gz(2*i+1)=gz(2*i+1)-zij*dphi/r12

      enddo

c     electric field

      do i=1,nH2+1
         Ex(i)=0.D0
         Ey(i)=0.D0
         Ez(i)=0.D0
      enddo

c     pair energies (+ gradient)

      ih=1
      
      do i=1,nH2
         do k=1,3
            xih=xhm-xH2(i,k)
            yih=yhm-yH2(i,k)
            zih=zhm-zH2(i,k)
            xih2=xih**2
            yih2=yih**2
            zih2=zih**2
            rih2=xih2+yih2+zih2
            rih=sqrt(rih2)
            rih3=rih2*rih

            ddamp=0.D0
            if (rih.lt.rcut) then
               damp=exp(-(1.D0-rcut/rih)**2)
               ddamp=2.D0*rcut*(1.D0-rcut/rih)
            else
               damp=1.D0
            endif
            
            ecoul=ecoul+damp*q(1)*qH2(k)/rih

            ecoul=ecoul+facrep/rih**12

               dum=damp*(1.D0+ddamp/rih)*q(1)*qH2(k)/rih3
               dum=dum+12.D0*facrep/rih**14

               if (k.eq.1) then
                  gx(2*i  )=gx(2*i  )+dum*xih
                  gy(2*i  )=gy(2*i  )+dum*yih
                  gz(2*i  )=gz(2*i  )+dum*zih

                  gx(1)=gx(1)-dum*xih
                  gy(1)=gy(1)-dum*yih
                  gz(1)=gz(1)-dum*zih
               elseif (k.eq.2) then
                  gx(2*i+1)=gx(2*i+1)+dum*xih
                  gy(2*i+1)=gy(2*i+1)+dum*yih
                  gz(2*i+1)=gz(2*i+1)+dum*zih

                  gx(1)=gx(1)-dum*xih
                  gy(1)=gy(1)-dum*yih
                  gz(1)=gz(1)-dum*zih
               else             ! k=3: affects both sites
                  gx(2*i  )=gx(2*i  )+dum*xih/2.D0
                  gy(2*i  )=gy(2*i  )+dum*yih/2.D0
                  gz(2*i  )=gz(2*i  )+dum*zih/2.D0
                  gx(2*i+1)=gx(2*i+1)+dum*xih/2.D0
                  gy(2*i+1)=gy(2*i+1)+dum*yih/2.D0
                  gz(2*i+1)=gz(2*i+1)+dum*zih/2.D0

                  gx(1)=gx(1)-dum*xih
                  gy(1)=gy(1)-dum*yih
                  gz(1)=gz(1)-dum*zih
               endif

c     field on ion

            Ex(1)=Ex(1)+damp*qH2(k)*xih/rih3
            Ey(1)=Ey(1)+damp*qH2(k)*yih/rih3
            Ez(1)=Ez(1)+damp*qH2(k)*zih/rih3

               dexx(1,3*i-3+k)=-damp*qH2(k)*(1.D0-3.D0*xih2/rih2)/rih3
               deyy(1,3*i-3+k)=-damp*qH2(k)*(1.D0-3.D0*yih2/rih2)/rih3
               dezz(1,3*i-3+k)=-damp*qH2(k)*(1.D0-3.D0*zih2/rih2)/rih3
               dexy(1,3*i-3+k)= damp*qH2(k)*3.D0*xih*yih/rih3/rih2
               dexz(1,3*i-3+k)= damp*qH2(k)*3.D0*xih*zih/rih3/rih2
               deyz(1,3*i-3+k)= damp*qH2(k)*3.D0*yih*zih/rih3/rih2
            
c     bare electric fields on H2 centers (site k=3)

            if (k.eq.3) then
               Ex(i+1)=Ex(i+1)-q(1)*xih/rih3
               Ey(i+1)=Ey(i+1)-q(1)*yih/rih3
               Ez(i+1)=Ez(i+1)-q(1)*zih/rih3

                  dexx(i+1,ih)=-q(1)*(1.D0-3.D0*xih2/rih2)/rih3
                  deyy(i+1,ih)=-q(1)*(1.D0-3.D0*yih2/rih2)/rih3
                  dezz(i+1,ih)=-q(1)*(1.D0-3.D0*zih2/rih2)/rih3
                  dexy(i+1,ih)= q(1)*3.D0*xih*yih/rih3/rih2
                  dexz(i+1,ih)= q(1)*3.D0*xih*zih/rih3/rih2
                  deyz(i+1,ih)= q(1)*3.D0*yih*zih/rih3/rih2
               
c     repulsion-dispersion (also if site k=3)

               erep=erep+A*exp(-b*rih)

               hcdmp=1.d0
               dhcdmp=0.d0
               if (rih.lt.rdmp) then
                  hcdmp=exp(-(rdmp/rih-1.d0)**2)
                     dhcdmp=2.d0*rdmp*(rdmp/rih-1.d0)*hcdmp/rih2
               endif
                  
               hcdisp=C6/rih2**3+C8/rih2**4
                  
               edisp=edisp-hcdisp*hcdmp
                  
                  dumrep=-b*A*(exp(-b*rih))/rih
                  dumdisp=(6.d0*C6/rih2**3+8.d0*C8/rih2**4)/rih2

                  dum=-dumdisp*hcdmp+dhcdmp*hcdisp/rih
                  
                  dum=(dum-dumrep)/2.D0
                  
                  gx(2*i  )=gx(2*i  )+xih*dum
                  gy(2*i  )=gy(2*i  )+yih*dum
                  gz(2*i  )=gz(2*i  )+zih*dum
                  
                  gx(2*i+1)=gx(2*i+1)+xih*dum
                  gy(2*i+1)=gy(2*i+1)+yih*dum
                  gz(2*i+1)=gz(2*i+1)+zih*dum

                  gx(1)=gx(1)-2D0*xih*dum
                  gy(1)=gy(1)-2D0*yih*dum
                  gz(1)=gz(1)-2D0*zih*dum
               
            endif

         enddo
      enddo

c     H2-H2 interaction

      do i=1,nH2
         do j=i+1,nH2
            do k=1,3
               do l=1,3
                  xik_jl=xH2(i,k)-xH2(j,l)
                  yik_jl=yH2(i,k)-yH2(j,l)
                  zik_jl=zH2(i,k)-zH2(j,l)
                  rik_jl_2=xik_jl**2+yik_jl**2+zik_jl**2
                  rik_jl=sqrt(rik_jl_2)
                  dumq=qH2(k)*qH2(l)/rik_jl

                  ecoul=ecoul+dumq

                  dumq=-dumq/rik_jl_2
                  
                     if (k.lt.3.and.l.lt.3) then
                        gx(2*i-1+k)=gx(2*i-1+k)+xik_jl*dumq
                        gy(2*i-1+k)=gy(2*i-1+k)+yik_jl*dumq
                        gz(2*i-1+k)=gz(2*i-1+k)+zik_jl*dumq

                        gx(2*j-1+l)=gx(2*j-1+l)-xik_jl*dumq
                        gy(2*j-1+l)=gy(2*j-1+l)-yik_jl*dumq
                        gz(2*j-1+l)=gz(2*j-1+l)-zik_jl*dumq
                     elseif (k.lt.3) then
                        gx(2*i-1+k)=gx(2*i-1+k)+xik_jl*dumq
                        gy(2*i-1+k)=gy(2*i-1+k)+yik_jl*dumq
                        gz(2*i-1+k)=gz(2*i-1+k)+zik_jl*dumq

                        gx(2*j  )=gx(2*j  )-xik_jl*dumq/2.D0
                        gy(2*j  )=gy(2*j  )-yik_jl*dumq/2.D0
                        gz(2*j  )=gz(2*j  )-zik_jl*dumq/2.D0

                        gx(2*j+1)=gx(2*j+1)-xik_jl*dumq/2.D0
                        gy(2*j+1)=gy(2*j+1)-yik_jl*dumq/2.D0
                        gz(2*j+1)=gz(2*j+1)-zik_jl*dumq/2.D0
                     elseif (l.lt.3) then
                        gx(2*j-1+l)=gx(2*j-1+l)-xik_jl*dumq
                        gy(2*j-1+l)=gy(2*j-1+l)-yik_jl*dumq
                        gz(2*j-1+l)=gz(2*j-1+l)-zik_jl*dumq

                        gx(2*i  )=gx(2*i  )+xik_jl*dumq/2.D0
                        gy(2*i  )=gy(2*i  )+yik_jl*dumq/2.D0
                        gz(2*i  )=gz(2*i  )+zik_jl*dumq/2.D0

                        gx(2*i+1)=gx(2*i+1)+xik_jl*dumq/2.D0
                        gy(2*i+1)=gy(2*i+1)+yik_jl*dumq/2.D0
                        gz(2*i+1)=gz(2*i+1)+zik_jl*dumq/2.D0
                     else
                        gx(2*i  )=gx(2*i  )+xik_jl*dumq/2.D0
                        gy(2*i  )=gy(2*i  )+yik_jl*dumq/2.D0
                        gz(2*i  )=gz(2*i  )+zik_jl*dumq/2.D0

                        gx(2*i+1)=gx(2*i+1)+xik_jl*dumq/2.D0
                        gy(2*i+1)=gy(2*i+1)+yik_jl*dumq/2.D0
                        gz(2*i+1)=gz(2*i+1)+zik_jl*dumq/2.D0

                        gx(2*j  )=gx(2*j  )-xik_jl*dumq/2.D0
                        gy(2*j  )=gy(2*j  )-yik_jl*dumq/2.D0
                        gz(2*j  )=gz(2*j  )-zik_jl*dumq/2.D0

                        gx(2*j+1)=gx(2*j+1)-xik_jl*dumq/2.D0
                        gy(2*j+1)=gy(2*j+1)-yik_jl*dumq/2.D0
                        gz(2*j+1)=gz(2*j+1)-zik_jl*dumq/2.D0
                     endif

                  if (k.eq.3.and.l.eq.3) then ! LJ between c.o.m.

                     rik_jl_6=(sigma2/rik_jl_2)**3
                     eH2=eH2+4.D0*epsilon*rik_jl_6*(rik_jl_6-1.D0)

                      dumLJ=-12.D0*epsilon*rik_jl_6*(2.D0*rik_jl_6-1.D0)
                        gx(2*i  )=gx(2*i  )+xik_jl*dumLJ/rik_jl_2
                        gy(2*i  )=gy(2*i  )+yik_jl*dumLJ/rik_jl_2
                        gz(2*i  )=gz(2*i  )+zik_jl*dumLJ/rik_jl_2

                        gx(2*i+1)=gx(2*i+1)+xik_jl*dumLJ/rik_jl_2
                        gy(2*i+1)=gy(2*i+1)+yik_jl*dumLJ/rik_jl_2
                        gz(2*i+1)=gz(2*i+1)+zik_jl*dumLJ/rik_jl_2

                        gx(2*j  )=gx(2*j  )-xik_jl*dumLJ/rik_jl_2
                        gy(2*j  )=gy(2*j  )-yik_jl*dumLJ/rik_jl_2
                        gz(2*j  )=gz(2*j  )-zik_jl*dumLJ/rik_jl_2

                        gx(2*j+1)=gx(2*j+1)-xik_jl*dumLJ/rik_jl_2
                        gy(2*j+1)=gy(2*j+1)-yik_jl*dumLJ/rik_jl_2
                        gz(2*j+1)=gz(2*j+1)-zik_jl*dumLJ/rik_jl_2

                  endif

               enddo
            enddo
         enddo
      enddo

c     polarization

      do i=1,nH2+1
         epol=epol-polar(i)*(Ex(i)**2+Ey(i)**2+Ez(i)**2)/2.D0
      enddo

c     contribution of charges on H2 to polarization on H-
         do i=1,nh2
            do k=1,3
               dumx=ex(1)*dexx(1,3*i-3+k)+ey(1)*dexy(1,3*i-3+k)+ez(1)
     $              *dexz(1,3*i-3+k)
               dumy=ex(1)*dexy(1,3*i-3+k)+ey(1)*deyy(1,3*i-3+k)+ez(1)
     $              *deyz(1,3*i-3+k)
               dumz=ex(1)*dexz(1,3*i-3+k)+ey(1)*deyz(1,3*i-3+k)+ez(1)
     $              *dezz(1,3*i-3+k)

               dumx=-dumx
               dumy=-dumy
               dumz=-dumz
               
               gxpol(1)=gxpol(1)-polar(1)*dumx
               gypol(1)=gypol(1)-polar(1)*dumy
               gzpol(1)=gzpol(1)-polar(1)*dumz

               if (k.eq.1) then
                  gxpol(2*i  )=gxpol(2*i  )+polar(1)*dumx
                  gypol(2*i  )=gypol(2*i  )+polar(1)*dumy
                  gzpol(2*i  )=gzpol(2*i  )+polar(1)*dumz
               elseif (k.eq.2) then
                  gxpol(2*i+1)=gxpol(2*i+1)+polar(1)*dumx
                  gypol(2*i+1)=gypol(2*i+1)+polar(1)*dumy
                  gzpol(2*i+1)=gzpol(2*i+1)+polar(1)*dumz
               else
                  gxpol(2*i  )=gxpol(2*i  )+polar(1)*dumx/2.D0
                  gypol(2*i  )=gypol(2*i  )+polar(1)*dumy/2.D0
                  gzpol(2*i  )=gzpol(2*i  )+polar(1)*dumz/2.D0
                  gxpol(2*i+1)=gxpol(2*i+1)+polar(1)*dumx/2.D0
                  gypol(2*i+1)=gypol(2*i+1)+polar(1)*dumy/2.D0
                  gzpol(2*i+1)=gzpol(2*i+1)+polar(1)*dumz/2.D0
               endif
               
            enddo
         enddo

c     contribution of H- to polarization of H2

         ih=1
         do i=1,nH2
            dumx=ex(i+1)*dexx(i+1,ih)+ey(i+1)*dexy(i+1,ih)+ez(i+1)
     $           *dexz(i+1,ih)
            dumy=ex(i+1)*dexy(i+1,ih)+ey(i+1)*deyy(i+1,ih)+ez(i+1)
     $           *deyz(i+1,ih)
            dumz=ex(i+1)*dexz(i+1,ih)+ey(i+1)*deyz(i+1,ih)+ez(i+1)
     $           *dezz(i+1,ih)
            
            dumx=-dumx
            dumy=-dumy
            dumz=-dumz

            gxpol(2*i  )=gxpol(2*i  )-polar(i+1)*dumx/2.D0
            gypol(2*i  )=gypol(2*i  )-polar(i+1)*dumy/2.D0
            gzpol(2*i  )=gzpol(2*i  )-polar(i+1)*dumz/2.D0
            
            gxpol(2*i+1)=gxpol(2*i+1)-polar(i+1)*dumx/2.D0
            gypol(2*i+1)=gypol(2*i+1)-polar(i+1)*dumy/2.D0
            gzpol(2*i+1)=gzpol(2*i+1)-polar(i+1)*dumz/2.D0
            
            gxpol(1)=gxpol(1)+polar(i+1)*dumx
            gypol(1)=gypol(1)+polar(i+1)*dumy
            gzpol(1)=gzpol(1)+polar(i+1)*dumz
         enddo

      epot=eintra+ecoul+erep+edisp+eH2+epol

      do i=1,Natoms
         dP(3*i-2)= (gx(i)+gxpol(i))
         dP(3*i-1)= (gy(i)+gypol(i))
         dP(3*i  )= (gz(i)+gzpol(i))
      enddo

      return
      end

