function[MO, BO]= SolvePiece(ev,  Pv, Tv, j, numSalti, x0, MOtemp, BOtemp, numStati, numStati_a, numStati_c, numStati_p, numStati_r, flag)
  MO=MOtemp;
  BO=BOtemp;
  if(j==0)
      E1=ev(:,1:numStati_a);
      P1=Pv(1:numStati,1:numStati);
      T=Tv(1:numStati,1:numStati_a);
      if(flag)
        newP=T'*P1*T-E1'*T'*P1*T*E1;
      else
        newP=Pv;
      end
      
      Rc=newP(1:numStati_c,1:numStati_c);
      Rpr=newP(1:numStati_c,numStati_c+1:end);
      xnc=x0(numStati_c+1:end,1);
    
      MO(1:numStati_c, 1:numStati_c)=MO(1:numStati_c, 1:numStati_c)+Rc;
    
      BO(1:numStati_c,1)=BO(1:numStati_c,1)-Rpr*xnc;
  else
        %extract the last e^At
        Ej=ev(:,(j-1)*numStati_a+1:j*numStati_a);
        % computed P2
        Pnext=Pv(1:numStati,numStati*(j)+1:numStati*(j+1));
        
        %extract Ec, Ep and Er pieces from Ej
        Ec=Ej(numStati_c+1:numStati_c+numStati_p, 1:numStati_c);
        Ep=Ej(numStati_c+1:numStati_c+numStati_p, numStati_c+1:numStati_c+numStati_p);
        Er=Ej(numStati_c+1:numStati_c+numStati_p, numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
        

        % computed T2
        Tnext=Tv(1:numStati,numStati_a*(j)+1:numStati_a*(j+1));        

        % compute the piece 2xcq P(q+q, cpp)epq xp(q-1)
        index_lambda=(numSalti+1)*(numStati_c);

        row_M=zeros(numStati_p, index_lambda);
        xnc=x0(numStati_c+1:end,1);
        chain_xp=xnc(1:numStati_p,1);
        xr_prev=xnc(numStati_p+1:end,1);
        
        for k=1:j
            % compute the piece 2xcq P(q+1, cpp)epq xp(q-1)
            Ek=ev(:,(k-1)*numStati_a+1:(k)*numStati_a);
            
            %extract Ec, Ep and Er pieces from Ek
            Eck=Ek(numStati_c+1:numStati_c+numStati_p, 1:numStati_c);
            Epk=Ek(numStati_c+1:numStati_c+numStati_p, numStati_c+1:numStati_c+numStati_p);
            Erk=Ek(numStati_c+1:numStati_c+numStati_p, numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
            row_M=Epk*row_M;
            row_M(:,(k-1)*numStati_c+1: k*numStati_c)=Eck;
            chain_xp=Ep*chain_xp+Er*xr_prev;
            
            xr_prev=Ek(numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r, numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r)*xr_prev;
        end

        
        % P2-eAt*P3*eAt
        % computes new dynamic if needed
        if(flag)
            newP=Tnext'*Pnext*Tnext;
            if((j+2)<=(numSalti+1))
                dynamics =ev(:,(k)*numStati_a+1:(k+1)*numStati_a);
                newP=newP-dynamics'*Tnext'*Pnext*Tnext*dynamics;
            end
        else
            newP=Pv;
        end
        
        
        chainP=newP;
        for k=j:-1:1
            %compute the piect 2xp(q-1)epq P(q+1, prqp) xrq
            
            % compute the piece 2xcq P(q+1, cpp)epq xp(q-1)
            Ek=ev(:,(k-1)*numStati_a+1:(k)*numStati_a);

            %extract Ec, Ep and Er pieces from Ek
            Eck=Ek(numStati_c+1:numStati_c+numStati_p, 1:numStati_c);
            Epk=Ek(numStati_c+1:numStati_c+numStati_p, numStati_c+1:numStati_c+numStati_p);
            BO(numStati_c*(k-1)+1:numStati_c*k,1)=BO(numStati_c*(k-1)+1: numStati_c*k,1)-Eck'*chainP(numStati_c+1:numStati_c+numStati_p,numStati_c+numStati_p+1:end)*xr_prev;
            chainP=Epk'*chainP(numStati_c+1:numStati_c+numStati_p,numStati_c+numStati_p+1:end);
        end
   
        %build this piece
        Rc=newP(1:numStati_c,1:numStati_c);
        Rpr=newP(1:numStati_c,numStati_c+1:end);
     
        MO(j*numStati_c+1:(j+1)*numStati_c,1:index_lambda)=MO(j*numStati_c+1:(j+1)*numStati_c,1:index_lambda)+Rpr(:,1:numStati_p)*row_M;
        MO(1:index_lambda,j*numStati_c+1:(j+1)*numStati_c)=MO(1:index_lambda,j*numStati_c+1:(j+1)*numStati_c)+ MO(j*numStati_c+1:(j+1)*numStati_c,1:index_lambda)';
        
        MO(j*numStati_c+1:(j+1)*numStati_c,j*numStati_c+1:(j+1)*numStati_c)=MO(j*numStati_c+1:(j+1)*numStati_c,j*numStati_c+1:(j+1)*numStati_c)+Rc;
        
        
        BO(j*numStati_c+1:(j+1)*numStati_c,1)=BO(j*numStati_c+1:(j+1)*numStati_c,1)-Rpr(:,1:numStati_p)*chain_xp-Rpr(:,numStati_p+1:end)*xr_prev;
    
        %build the matrix for the recursive
        Px=[Ec, Ep, Er];
        Prec=Px'*newP(numStati_c+1:numStati_c+numStati_p,numStati_c+1:numStati_c+numStati_p)*Px;
        %go recursively on this
        [MO, BO]=SolvePiece(ev,  Prec, Tv, j-1, numSalti, x0, MO, BO, numStati, numStati_a, numStati_c, numStati_p, numStati_r,0);
    
  end
end