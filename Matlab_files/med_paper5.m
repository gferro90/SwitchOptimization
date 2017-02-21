clc 
clear all
close all

% Aggiunge il FF al med_paper1 e ottimizza [CI_c CI_w] (ancora senza
% formula quindi non va bene)


%% DEFINIZIONE TEMPO DI SIMULAZIONE

% Definizione di un predefinito numero di campioni spalmati sull'intervallo
% di tempo
numeroCampioni=20000;

% Definizione del tempo finale di simulazione
tempoFinale=2;

% Definizione del delta t
delta_t=tempoFinale/numeroCampioni;

numeroCampioni=numeroCampioni+1;

% Definizione del vettore dei tempi
t=[0:delta_t:tempoFinale];

%% INPUT DATI
% Input da questo file.
code=0;

if(code==1)
    %inserire il numero di salti:
    numControllori=input('Inserire il numero desiderato di Controllori: ');
    numSalti=numControllori-1;
    
    %inserire il riferimento (quello a regime permanente):
    numRifStatico=input('Inserire il numeratore del riferimento statico: ');
    denRifStatico=input('Inserire il denominatore del riferimento statico: ');
    rifStaticoTempo=input('Inserire l espressione del riferimento statico nel tempo: ');
    
    %inserire il riferimento dinamico:
    numRifDinamico=input('Inserire il numeratore del riferimento dinamico: ');
    denRifDinamico=input('Inserire il denominatore del riferimento dinamico: ');
    rifDinamicoTempo=input('Inserire l espressione del riferimento dinamico nel tempo: ');

    
    %inserire il processo:
    numProcesso=input('Inserire il numeratore del processo: ');
    
    %inserire condizioni iniziali del processo:
    condInizialiProcesso=input('Inserire le condizioni iniziali del processo: ');
    
    for i=1:numControllori
        stringaControllerNum=strcat('Inserire il numeratore del controllore ',num2str(i),': ');
        controllerNum=input(stringaControllerNum);
        stringaControllerDen=strcat('Inserire il denominatore del controllore ',num2str(i),': ');
        controllerDen=input(stringaControllerDen);
        sistema(i).fdtControllore=tf(controllerNum, controllerDen);
    end
    
    guadagnoErrore=input('Inserire guadagno relativo all errore: ');
    guadagnoDerivataErrore=input('Inserire guadagno relativo alla derivata dell errore');
    
else
    numSalti=1;
    numControllori=numSalti+1;

    esponente=20;
    
     numRifDinamico=esponente;
     denRifDinamico=[1 esponente 0];
%     numRifDinamico=1;
%     denRifDinamico=[1 0];
    rifDinamicoTempo=1-exp(-esponente*t);
    
    wn=2*pi/0.3;
    delta=0.8;
    % P=k/(s^2+as+b)
    k=36*(wn^2);
    a=2*delta*wn;
    b=wn^2;
    
    numProcesso=[k];
    denProcesso=[1 a b];
    
    CI_Processo=[0; 0];
    CI_Controllore=0;
    CI_cc=[CI_Controllore; CI_Processo];

    Kp=0.1;
    Ki=[0.2 3];
   % Ki=[0.2 3 0.2];

    sistema(1).fdtControllore=tf([Kp Ki(1)], [1 0]);
    sistema(2).fdtControllore=tf([Kp Ki(2)], [1 0]);
   % sistema(3).fdtControllore=tf([Kp Ki(3)], [1 0]);
   
   % segnale di feedforward
   den_w=[1 100];
   num_w=[1];
 
    % guadagno associato all'errore nel funzionale di costo
    p1=10;
    % guadagno associato alla derivata dell'errore nel funzionale di costo
    p2=1e-2;
end

%% COSTRUZIONE FUNZIONI DI TRASFERIMENTO

%fdt del riferimento dinamico
fdtRifDinamico=tf(numRifDinamico, denRifDinamico)

%fdt del processo
fdtProcesso=tf(numProcesso, denProcesso)

%fdt del FF
fdtFF=tf(num_w, den_w);


%fdt dei sistemi a ciclo chiuso 
for i=1:numControllori
    %il sistema a ciclo chiuso
    sistema(i).fdtCicloChiuso=feedback(sistema(i).fdtControllore*fdtProcesso,1)
    sistema(i).fdtCicloChiusoFF=feedback(fdtProcesso,sistema(i).fdtControllore)
end


%% COSTRUZIONE MATRICI DI STATO
% phi = rd - rs
% COSTRUZIONE A_phi
   
% Passare al sistema in forma minima per semplificare la funzione di
% trasferimento
sistema_rd=ss(fdtRifDinamico,'minimal');
   
% Ottenere numeratore e denominatore della fdt
[num_rd,den_rd]=tfdata(tf(sistema_rd),'v');
    
% Calcolo del numero di stati
numStati_r=length(den_rd)-1;

A_rd=zeros(numStati_r, numStati_r);

% Calcolo della matrice A del sistema
if(numStati_r>1)
    A_rd(1:numStati_r-1,2:numStati_r)=eye(numStati_r-1);
    A_rd(1:numStati_r-1,1)=zeros(numStati_r-1,1);
    A_rd(numStati_r,:)=-den_rd(numStati_r+1:-1:2);
else
    if (numStati_r==1)
        A_rd=-den_rd(2);
    end
end
 % Costruzione del sistema processo
[num_p,den_p]=tfdata(fdtProcesso,'v');
numStati_p=length(den_p)-1;
if(numStati_p>1)
    % 1) Riempire le prime n-1 righe dalla seconda colonna in poi con I shiftato a destra di 1.
    A_p(1:numStati_p-1,2:numStati_p)=eye(numStati_p-1);
    % 2) Riempire la prima colonna fino a n-1 righe di zeri.
    A_p(1:numStati_p-1,1)=zeros(numStati_p-1,1);
    % 3) Riempire l'ultima riga con i coefficienti del denominatore
    % cambiati di segno.
    A_p(numStati_p,:)=-den_p(numStati_p+1:-1:2);
else
    % In questo caso A è scalare.
    A_p=-den_p(2);
end

% Costruire B in forma canonica di controllore:
% 1) Riempire le prime n-1 righe di zeri.
B_p=zeros(numStati_p,1);
% 2) Settare a 1 l'ultima riga.
B_p(numStati_p,1)=1;

D_p=num_p(1);

% Costruire C in forma canonica di controllore:
C_p=zeros(1,numStati_p);
% 1) Riempire la prima riga di C con i coefficienti del numeratore dal
% grado + basso al + alto.
C_p(1:numStati_p)=num_p(numStati_p+ 1:-1:2)-D_p*den_p(numStati_p+ 1:-1:2);

% Normalizzo la matrice C mettendo i suoi coefficienti nella dinamica
mul=eye(numStati_p);
for j=1:numStati_p
    if(C_p(j)~=0)
        mul(j,j)=C_p(j);
        C_p(j)=1;
    end
end
A_p=mul*A_p*(mul^-1);
B_p=mul*B_p;

% Passare al sistema in forma minima per semplificare la funzione di
% trasferimento
sistema_w=ss(fdtFF,'minimal');
   
% Ottenere numeratore e denominatore della fdt
[num_w,den_w]=tfdata(tf(sistema_w),'v');
%Calcolo del numero di stati
numStati_w=length(den_w)-1;

A_w=zeros(numStati_w, numStati_w);

% Calcolo della matrice A del sistema
if(numStati_w>1)
    A_w(1:numStati_w-1,2:numStati_w)=eye(numStati_w-1);
    A_w(1:numStati_w-1,1)=zeros(numStati_w-1,1);
    A_w(numStati_w,:)=-den_w(numStati_w+1:-1:2);
else
    if (numStati_w==1)
        A_w=-den_w(2);
    end
end
 

Ifull=eye(numStati_r);
Ir=Ifull(1,:);


Ifull=eye(numStati_w);
Iw=Ifull(1,:);

B_w=B_p;
D_w=D_p;
for i=1:numControllori
     
    % Ricava numeratore e denominatore della fdt a ciclo chiuso. 
    % I vettori sono di uguale lunghezza e vanno dal grado maggiore a
    % quello minore.
    [num_cc, den_cc]=tfdata(sistema(i).fdtCicloChiuso,'v')
    
    % Ricava il numero di stati a ciclo chiuso uguale al grado del
    % denominatore della fdt.
    numStati_cc=length(den_cc)-1;
    
    % Costruire Acc considerando come stato [xc xp]
    % Costruzione del sistema controllore 
    [num_c,den_c]=tfdata(sistema(i).fdtControllore,'v');

    numStati_c=length(den_c)-1;
    A_c=zeros(numStati_c, numStati_c);
    if(numStati_c>1)
        % 1) Riempire le prime n-1 righe dalla seconda colonna in poi con I shiftato a destra di 1.
        A_c(1:numStati_c-1,2:numStati_c)=eye(numStati_c-1);
        % 2) Riempire la prima colonna fino a n-1 righe di zeri.
        A_c(1:numStati_c-1,1)=zeros(numStati_c-1,1);
        % 3) Riempire l'ultima riga con i coefficienti del denominatore
        % cambiati di segno.
        A_c(numStati_c,:)=-den_c(numStati_c+1:-1:2);
    else
        % In questo caso A è scalare.
        A_c=-den_c(2);
    end

    
    % Costruire B in forma canonica di controllore:
    % 1) Riempire le prime n-1 righe di zeri.
    B_c=zeros(numStati_c,1);
    % 2) Settare a 1 l'ultima riga.
    B_c(numStati_c,1)=1;
    
    % Costruire C in forma canonica di controllore:
    C_c=zeros(1,numStati_c);
    
    
    % Se != 0 D è il coefficiente del termine di grado massimo
    D_c=num_c(1);
    
    % 1) Riempire la prima riga di C con i coefficienti del numeratore dal
    % grado + basso al + alto.
    C_c(1:numStati_c)=num_c(numStati_c+ 1:-1:2)- D_c*den_c(numStati_c+ 1:-1:2);
    
           
    % Normalizzo la matrice C mettendo i suoi coefficienti nella dinamica
    mul=eye(numStati_c);
    for j=1:numStati_c
        if(C_c(j)~=0)
            mul(j,j)=C_c(j);
            C_c(j)=1;
        end
    end
    A_c=mul*A_c*(mul^-1);
    B_c=mul*B_c;
    
    
    sistema(i).A_c=A_c;
    sistema(i).B_c=B_c;
    sistema(i).C_c=C_c;
    sistema(i).D_c=D_c;
    
    L=1+D_p*D_c;
    
   
    
    % Costruzione sistema a ciclo chiuso
    sistema(i).A_cc=[A_c-B_c*(L^-1)*D_p*C_c               -B_c*(L^-1)*C_p; 
                     B_p*C_c-B_p*D_c*(L^-1)*D_p*C_c       A_p-B_p*D_c*(L^-1)*C_p];
    sistema(i).B_cc=[B_c-B_c*(L^-1)*D_p*D_c; 
                     B_p*D_c-B_p*D_c*(L^-1)*D_p*D_c];
                 
    sistema(i).B_cw=[-B_c*(L^-1)*D_p; 
                     -B_p*D_c*(L^-1)*D_p+B_w];
    sistema(i).C_cc=[(L^-1)*D_p*C_c (L^-1)*C_p];
    sistema(i).D_cc=(L^-1)*D_p*D_c;
    sistema(i).D_cw=D_w;
    
    disp(eig(sistema(i).A_cc));
    
    sistema(i).A_a=[sistema(i).A_cc,               sistema(i).B_cc*Ir,        sistema(i).B_cw*Iw;
        zeros(numStati_r,numStati_cc),                   A_rd,                zeros(numStati_r, numStati_w);
        zeros(numStati_w, numStati_cc),      zeros(numStati_w, numStati_r),                A_w];
    
    
    if(i==1)
        % La risposta libera del sistema a ciclo chiuso nel dominio di Laplace
        fdtRispostaLiberaCicloChiuso=tf(ss(sistema(i).A_cc, CI_cc, sistema(i).C_cc, sistema(i).D_cc));
    else
        fdtRispostaLiberaCicloChiuso=0;
    end
    
  
    
    %sistema(i).T_normal=[-sistema(i).O (eye(numStati)-sistema(i).M1)*Nr ];
    % L'errore in evoluzione libera
    sistema(i).fdt_ed=fdtRifDinamico-(sistema(i).fdtCicloChiuso*fdtRifDinamico);
    %-fdtRispostaLiberaCicloChiuso;

    % Passare al sistema in forma minima per semplificare la funzione di
    % trasferimento
    sistema_ed=ss(sistema(i).fdt_ed,'minimal');

    % Ottenere numeratore e denominatore della funzione in Laplace di y-rs
    [num_ed, den_ed]=tfdata(tf(sistema_ed), 'v');
    
    % Calcolare la dimensione della matrice Ai
    numStati_ed=length(den_ed)-1;
    
    % Salvare i dati del primo sistema che serviranno poi per il calcolo
    % delle condizioni iniziali
    if(i==1)
        num_ed_i=num_ed;
        den_ed_i=den_ed;
        numStati_ed_i=numStati_ed;
    end
    
    % Costruzione della matrice Ai del sistema in evoluzione libera
    % ded=Ai*ed con ed=rd-y
    A_ed_temp=zeros(numStati_ed, numStati_ed);
    
    if(numStati_ed>1)
        A_ed_temp(1:numStati_ed-1, 2:numStati_ed)=eye(numStati_ed-1);
        A_ed_temp(1:numStati_ed-1, 1)=zeros(numStati_ed-1,1);
        A_ed_temp(numStati_ed,:)=-den_ed(numStati_ed+1:-1:2);
    else
        A_ed_temp=-den_ed(2);
    end
    
    numStati_v=numStati_ed;
      % calcolo matrice x il successivo cambio di coordinate a [e e' e'' ...]
    sistema(i).M1=zeros(numStati_v, numStati_v);
    sistema(i).M2=zeros(numStati_v, numStati_v);
    sistema(i).O=zeros(numStati_v, numStati_cc);
    
    
    % chain_1 is A^n    
    chain_1=eye(numStati_cc);
    % chain_2 is [A^(n-1)B | A^(n-2)B |...| B]
    chain_2=zeros(numStati_cc,numStati_v);
    chain_3=zeros(numStati_cc,numStati_v);
    for j=1:numStati_v
        sistema(i).O(j,1:numStati_cc)=sistema(i).C_cc*(chain_1);
        sistema(i).M1(j,:)=sistema(i).C_cc*(chain_2);
        sistema(i).M1(j,j)=sistema(i).D_cc;
        
        sistema(i).M2(j,:)=sistema(i).C_cc*(chain_3);
        sistema(i).M2(j,j)=sistema(i).D_cw;
        
        chain_1=sistema(i).A_cc*chain_1;
        chain_2=sistema(i).A_cc*chain_2;
        chain_2(:,j)=sistema(i).B_cc;
        
        chain_3=sistema(i).A_cc*chain_3;
        chain_3(:,j)=sistema(i).B_cw;
    end
    %sistema(i).M=eye(numStati, numStati)-sistema(i).M;
    % Mettere I per lo stato dell' "errore" dinamico ex_rd = rs-rd
    % alla fine facendo [x_r | 0] - M [x_cc | x_r | ex_rd] = [e de ... ex_rd]
    
        Nr=zeros(numStati_v, numStati_r);
    for kk=1:numStati_v
       Nr(kk,:)=Ir*(A_rd^(kk-1));
    end
    
    
    Nw=zeros(numStati_v, numStati_w);
    for kk=1:numStati_v
       Nw(kk,:)=Iw*(A_w^(kk-1));
    end
    
    % Calcolo della matrice T per questo sistema
    sistema(i).T_r=[-sistema(i).O, (eye(numStati_v)-sistema(i).M1)*Nr];
    sistema(i).T=[-sistema(i).O, (eye(numStati_v)-sistema(i).M1)*Nr, -sistema(i).M2*Nw;
                  zeros(numStati_w,numStati_cc), zeros(numStati_w, numStati_r) eye(numStati_w)];
    
    
   % B_a=[sistema(i).B_cw; zeros(numStati_r,1) ];
    B_ed=A_ed_temp*sistema(i).M2*Nw+sistema(i).T_r*[sistema(i).B_cw; zeros(numStati_r, numStati_w)]-sistema(i).M2*Nw*A_w;
    
    sistema(i).A_ed=[A_ed_temp,                        B_ed;
                     zeros(numStati_w,numStati_ed)     A_w];
                 
    numStati_ed=size(sistema(i).A_ed,1);

    numStati=numStati_ed;
       
    % numero stati iniziali
    numStati_i=numStati_ed_i;

    % Definizione del parametro epsilon per pesare gli stati non necessari
    epsilon=0;
    
    % Definizione della matrice Q usata nel funzionale di costo
    %sistema(i).Q=zeros(numStati, numStati);
    sistema(i).Q=eye(numStati, numStati)*epsilon;

    
    % Inserire i guadagni che pesano l'errore dinamico e la sua derivata
    sistema(i).Q(1,1)=p1;
    sistema(i).Q(2,2)=p2;
    
    % Calcolo della matrice P utilizzata nella legge di controllo switch
    sistema(i).P=lyap((sistema(i).A_ed)',sistema(i).Q);
    sistema(i).Px=sistema(i).T'*sistema(i).P*sistema(i).T;
    Kc_c=sistema(i).Px(1:numStati_c,1:numStati_c);
    Kp_c=sistema(i).Px(1:numStati_c,numStati_c+1:numStati_c+numStati_p);
    Kr_c=sistema(i).Px(1:numStati_c,numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
    Kw_c=sistema(i).Px(1:numStati_c,end-numStati_w+1:end);
    
    Kc_w=Kw_c';
    Kp_w=sistema(i).Px(end-numStati_w+1:end,numStati_c+1:numStati_c+numStati_p);
    Kr_w=sistema(i).Px(end-numStati_w+1:end,numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
    Kw_w=sistema(i).Px(end-numStati_w+1:end,end-numStati_w+1:end);
    
    
    sistema(i).Kcw=[Kc_c, Kw_c;
        Kc_w,  Kw_w ];
    sistema(i).Kp=[Kp_c;
        Kp_w];
    sistema(i).Kr=[Kr_c;
        Kr_w];
    
%     if(i>1)
%         % Calcolo matrici per cambio di coordinate da xcc a xe
%         sistema(i-1).F=sistema(i-1).O'*sistema(i).P*sistema(i-1).O;
%         sistema(i-1).G=sistema(i-1).O'*sistema(i).P*sistema(i-1).M;
%     end
end


%CI_w=0.02
%% CALCOLO DELLE CONDIZIONI INIZIALI DEL SISTEMA


% Calcolo della matrice M
M_ParametriCI_rd=zeros(numStati_r);

% Costruzione del vettore dei termini noti [a0 a1 ... a_n-1]
R_TerminiNotiCI_rd=(num_rd(numStati_r+1:-1:2))';

% Construzione della matrice dei parametri
for i=1:numStati_r
    % i coefficienti del denominatore vanno da numStati a 1
    M_ParametriCI_rd(i,1:numStati_r+1-i)=den_rd(numStati_r+1-i:-1:1)';
end
CI_rd=zeros(numStati_r,1);
% Calcolo delle condizioni iniziali come vettore colonna
CI_rd(1:numStati_r,1)=(M_ParametriCI_rd^-1)*R_TerminiNotiCI_rd
% for kk=numStati_r+1:numStati
%     CI_rd(kk,1)=Nr(kk-1, :)*CI_rd;
% end
      
% CALCOLO CONDIZIONI INIZIALI DI w (FF)
% Calcolo della matrice M
M_ParametriCI_w=zeros(numStati_w);

% Costruzione del vettore dei termini noti [a0 a1 ... a_n-1]
R_TerminiNotiCI_w=(num_w(numStati_w+1:-1:2))';

% Construzione della matrice dei parametri
for i=1:numStati_w
    % i coefficienti del denominatore vanno da numStati a 1
    M_ParametriCI_w(i,1:numStati_w+1-i)=den_w(numStati_w+1-i:-1:1)';
end
CI_w=zeros(numStati_w,1);
% Calcolo delle condizioni iniziali come vettore colonna
CI_w(1:numStati_w,1)=(M_ParametriCI_w^-1)*R_TerminiNotiCI_w


% CALCOLO CONDIZIONI INIZIALI DEL SISTEMA ed=rd-y

% Calcolo della matrice M
M_ParametriCI_ed=zeros(numStati_ed_i);

% Costruzione del vettore dei termini noti [a0 a1 ... a_n-1]
R_TerminiNotiCI_ed=(num_ed_i(numStati_ed_i+1:-1:2))';

% Construzione della matrice dei parametri
for i=1:numStati_ed_i
    % i coefficienti del denominatore vanno da numStati a 1
    M_ParametriCI_ed(i,1:numStati_ed_i+1-i)=den_ed_i(numStati_ed_i+1-i:-1:1)';
end
CI_w=0;
% Calcolo delle condizioni iniziali come vettore colonna
CI_ed=[(M_ParametriCI_ed^-1)*R_TerminiNotiCI_ed; CI_w]



%% CALCOLO TEMPI DI SWITCHING

% Definizione di un vettore per contenere gli istanti di switching grande
% al massimo quanto il numero di controllori meno uno
istantiSwitch=zeros(1, numSalti);

% Definizione di un vettore per contenere il funzionale di costo e la sua
% derivata

% costo con ottimizzazione sia del tempo di switch che dello stato di switch
J_optim=zeros(1,numeroCampioni);

% costo con solo ottimizzazione del tempo di switch
J_normal=zeros(1,numeroCampioni);

% Definizione del funzionale di costo iniziale
J_CostoIniziale=(CI_ed')*sistema(1).P*CI_ed;

% Condizione iniziale dello stato del riferimento
x_rd=CI_rd;

% Condizione iniziale dello stato del FF
x_w=CI_w;

% Condizione iniziale dopo ogni istante di switch. Inizialmente è la
% condizione iniziale del sistema precedentemente calcolata
CI_Salto=CI_ed;

% Condizione iniziale dello stato a ciclo chiuso
x_cc_optim=CI_cc;
x_cc_normal=CI_cc;

% Variabile per salvare l'istante dell'ultimo switch avvenuto per
% l'evoluzione dello stato del tempo
istanteUltimoSwitch=0;

% indice del controllore in uso pari al numero di switch +1
index_optim=1;
index_normal=1;
Jindex_optim=1;
Jindex_normal=1;

% flag per definire quando gli switch ammessi sono terminati
done_optim=0;
done_normal=0;

% inizializzazione derivata del funzionale di costo
dJ_optim=-1;
dJ_normal=-1;
dJ_no_optim=-1;

% storing dei dati per i plot
x_cc_optim_plot=zeros(1,numeroCampioni);
x_cc_normal_plot=zeros(1,numeroCampioni);

% evoluzione dei sistemi autonomi
d_x_cc_s=zeros(numStati_cc,numControllori);
x_cc_s=zeros(numStati_cc,numControllori);
for nSist=1:numControllori
    x_cc_s(:,nSist)=CI_cc;
end
x_cc_system_plot=zeros(numControllori, numeroCampioni);
x_rd_plot=zeros(numStati_r, numeroCampioni);
x_ed=zeros(numStati, 1);
x_ed_plot=zeros(numStati, numeroCampioni);
x_ed(:,1)=CI_ed;
x_w_plot=zeros(numStati_w, numeroCampioni);



for i=1:numeroCampioni
   
    % Calcolo delle matrici necessarie al calcolo dello stato xc ottimo.
%     Fa=sistema(Jindex_optim).F(1:numStati_c, 1:numStati_c);
%     Fb=sistema(Jindex_optim).F(1:numStati_c, numStati_c+1:numStati_cc);
%     Ga=sistema(Jindex_optim).G(1:numStati_c, 1:numStati);
     
    % estrazione dello stato del processo
    x_p_ts=x_cc_optim(numStati_c+1:numStati_cc,1);
    
    % calcolo dello stato ottimo del controllore  
    x_cw_ts=-(sistema(Jindex_optim+1).Kcw^-1)*(sistema(Jindex_optim+1).Kp*x_p_ts+sistema(Jindex_optim+1).Kr*x_rd);
    %-(Fa^(-1))*(Fb*x_p_ts-Ga*x_rd);
    
    % stato a ciclo chiuso del prossimo sistema
    x_cc_next=x_cc_optim;
    x_cc_next(1:numStati_c,1)=x_cw_ts(1:numStati_c);
    x_w_next=x_cw_ts(numStati_c+1:end);

    
%    S=[-sistema(Jindex_optim).O  sistema(Jindex_optim).M];
    xa_optim=[x_cc_optim; x_rd; x_w];
    xa_next=[x_cc_next; x_rd; x_w_next];
    xa_normal=[x_cc_normal; x_rd; zeros(numStati_w,1)];
    
    % stato aggregato [e_cc, x_(r-rd)] per entrambi i sistemi
    % stato finale x_e usato in J per entrambi i sistemi
    x_e_next=sistema(Jindex_optim+1).T*xa_next;
    x_e=sistema(Jindex_optim).T*xa_optim;
    x_e_normal=sistema(Jindex_normal).T*xa_normal;
    
    K_opt_cw=[-(sistema(Jindex_optim+1).Kcw^-1)*sistema(Jindex_optim+1).Kp*B_p*sistema(Jindex_optim).C_c,   -(sistema(Jindex_optim+1).Kcw^-1)*sistema(Jindex_optim+1).Kp*(A_p-B_p*sistema(Jindex_optim).D_c*C_p),   -(sistema(Jindex_optim+1).Kcw^-1)*(sistema(Jindex_optim+1).Kp*B_p*sistema(Jindex_optim).D_c*Ir+sistema(Jindex_optim+1).Kr*A_rd)    -(sistema(Jindex_optim+1).Kcw^-1)*sistema(Jindex_optim+1).Kp*B_w*Iw];
    
    
    % manca Im
    K_opt=[K_opt_cw(1:numStati_c,:) ;
    B_p*sistema(Jindex_optim).C_c,    A_p-B_p*sistema(Jindex_optim).D_c*C_p,    B_p*sistema(Jindex_optim).D_c*Ir, B_w*Iw;       
    zeros(numStati_r,numStati_c),    zeros(numStati_r,numStati_p),    A_rd,   zeros(numStati_r, numStati_w);
    K_opt_cw(numStati_c+1:end,:) ];

%     K=[-Fa^-1*Fb*B_p*sistema(Jindex_optim).C_c,  -Fa^-1*Fb*(A_p-B_p*sistema(Jindex_optim).D_c*C_p),  Fa^-1*(Ga*Nr-Fb*B_p*sistema(Jindex_optim).D_c*Ir);
%         B_p*sistema(Jindex_optim).C_c,                  A_p-B_p*sistema(Jindex_optim).D_c*C_p ,               B_p*sistema(Jindex_optim).D_c*Ir; 
%           zeros(numStati,numStati_c)                         zeros(numStati,numStati_p),                                Nr                          ];

    
    K_normal =[sistema(Jindex_normal).A_cc,             sistema(Jindex_normal).B_cc*Ir,                 sistema(Jindex_normal).B_cw*Iw;
            zeros(numStati_r, numStati_c+numStati_p),            A_rd,                                  zeros(numStati_r, numStati_w)    
            zeros(numStati_w, numStati_c+numStati_p)     zeros(numStati_w, numStati_r)                                   A_w];
    
    % Calcolo del funzionale di costo
    J_optim(i)=(x_e_next')*sistema(Jindex_optim+1).P*x_e_next-(x_e')*sistema(Jindex_optim).P*x_e+J_CostoIniziale;
    J_normal(i)=(x_e_normal')*sistema(Jindex_normal+1).P*x_e_normal-(x_e_normal')*sistema(Jindex_normal).P*x_e_normal+J_CostoIniziale;
    
    % CONSIDERARE QUI C_cc*x_cc
    x_cc_optim_plot(i)=sistema(index_optim).C_cc*x_cc_optim;
    x_cc_normal_plot(i)=sistema(index_normal).C_cc*x_cc_normal;
    x_rd_plot(:,i)=x_rd;
    x_ed_plot(:,i)=x_ed;
    x_w_plot(:,i)=x_w;
    
    % Calcolo della derivata del funzionale di costo. 
    if(i>1)
        dJ_optim=(J_optim(i)-J_optim(i-1))/delta_t;
        dJ_optim=xa_optim'*K_opt'*sistema(Jindex_optim+1).T'*sistema(Jindex_optim+1).P*sistema(Jindex_optim+1).T*xa_next+...
            xa_next'*sistema(Jindex_optim+1).T'*sistema(Jindex_optim+1).P*sistema(Jindex_optim+1).T*K_opt*xa_optim+...
            xa_optim'*sistema(Jindex_optim).T'*sistema(Jindex_optim).Q*sistema(Jindex_optim).T*xa_optim;
        dJ_normal=(J_normal(i)-J_normal(i-1))/delta_t;
        dJ_normal=xa_normal'*(K_normal'*sistema(Jindex_normal).T'*sistema(Jindex_normal+1).P*sistema(Jindex_normal).T + ...
            sistema(Jindex_normal).T'*sistema(Jindex_normal+1).P*sistema(Jindex_normal).T*K_normal+...
            sistema(Jindex_normal).T'*sistema(Jindex_normal).Q*sistema(Jindex_normal).T)*xa_normal;
    end
        
    if((dJ_optim>0) && (done_optim==0))
        if(J_optim(i)<J_CostoIniziale)
            % Salvataggio dell' istante di switch
            istantiSwitch(index_optim)=t(i);
            reset=x_cc_next
            if(index_optim>=numSalti)
                done_optim=1;
            else
                Jindex_optim=Jindex_optim+1;
            end
            x_cc_optim=x_cc_next;
            x_w=x_w_next;
            index_optim=index_optim+1;
            %x_ed_plot(:,i)=x_e_next;

            x_ed=x_e_next
           % J_CostoIniziale=(x_e_next')*sistema(Jindex_optim).P*x_e_next;
        end
    end
    
    
    
    if((dJ_normal>0) && (done_normal==0))
        % Salvataggio dell' istante di switch
        istantiSwitch_normal(index_normal)=t(i);
        
        if(index_normal>=numSalti)
            done_normal=1;
        else
            Jindex_normal=Jindex_normal+1;
        end
        
        index_normal=index_normal+1;
        
    end
       
   
    % Evoluzione dello stato a ciclo chiuso
    d_x_cc_optim=sistema(index_optim).A_cc*x_cc_optim+sistema(index_optim).B_cc*x_rd(1,1)+ sistema(index_optim).B_cw*x_w(1,1);
    x_cc_optim=x_cc_optim+d_x_cc_optim*delta_t;
            
    % evolutione dello stato senza ottimizzazione con switch al tempo
    % ottimo
    d_x_cc_normal=sistema(index_normal).A_cc*x_cc_normal+sistema(index_normal).B_cc*x_rd(1,1);
    x_cc_normal=x_cc_normal+d_x_cc_normal*delta_t;
    
    % Evoluzione dei sistemi 
    for nSist=1:numControllori
        d_x_cc_s(:,nSist)=sistema(nSist).A_cc* x_cc_s(:,nSist)+sistema(nSist).B_cc*x_rd(1,1);
        x_cc_s(:,nSist)=x_cc_s(:,nSist)+d_x_cc_s(:,nSist)*delta_t;
        x_cc_system_plot(nSist,i)=sistema(nSist).C_cc*x_cc_s(:,nSist);
    end
    
    % Evoluzione di rd
    d_x_rd=A_rd*x_rd;
    x_rd=x_rd+d_x_rd*delta_t;
    
    d_x_ed=sistema(index_optim).A_ed*x_ed;
    x_ed=x_ed+d_x_ed*delta_t;
    
    d_x_w=A_w*x_w;
    x_w=x_w+d_x_w*delta_t;
end

% Visualizza gli istanti di switch
istantiSwitch
istantiSwitch_normal

  
    
%% PLOTS

Plot1=figure(2);
ax(1)=subplot(2,1,1);
plot(t,x_cc_optim_plot,t,x_cc_normal_plot,t,x_rd_plot(1,:), t, x_cc_system_plot(1,:), t, x_cc_system_plot(2,:), t, x_ed_plot(1,:), t, x_w_plot(1,:));
grid on
h=legend('$C_{reset}$','$C_{switch}$', '$\overline{r}$', '$C_1$', '$C_2$', 'ed', 'w'); 
set(h,'Interpreter','latex')
plotParam

ax(2)=subplot(2,1,2);
plot(t,J_optim,'b',t,J_normal, 'r');
legend('Funzionale di Costo Reset', 'Funzionale di Costo Normal');   
grid on
plotParam   

linkaxes(ax,'x');