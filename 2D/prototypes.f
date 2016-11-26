



! p = degre
! problemType = euler, navier-stokes,...
! tend = temps de fin
! CI = champ d'inconnues initial
! CL = conditions aux limites
! mesh = maillage

function DFR(p,problemType,tend,CI,CL,mesh)
	! Fonction principale
	
	! SP : position des SP en xi (meme en eta) 							-> vecteur (P+1)
	! shapeFunctions : fonctions de forme aux SP et aux interfaces		-> tableau 3D 4x(P+3)x(P+3)
	! dShapeFunctions : derivees des fonctions de forme					-> 2 tableaux 3D 4x(P+3)x(P+3)
	
	! (Voir noms des grandeurs de maillage)
	! meshType: detecte le type du maillage (quadrangles,triangles,mixte)
	! Pour les conditions aux limites, voir avec les fonctions uInterface et fluxInterface ; voir aussi la structure de donnees du maillage
	
	
	! lagrangeP : polynomes de Lagrange de degre P							-> matrice (P+1)x(P+3) 
	! lagrangeP2 : polynomes de Lagrange de degre P+2						-> matrice (P+3)x(P+3)
	! dlagrangeP2 : derivees des polynomes de Lagrange de degre P+2			-> matrice (P+3)x(P+3)
	
	! (Jacobiennes ?)
	
	
	! (Initialisations des vecteurs contenant les inconnues)
	! tau, parametre de penalisation
	! alpha, coefficients de RK
end DFR




! uSPCell : inconnue aux SP d'une cellule
! lagrangeP : polynomes de Lagrangede  degre P
function extrapolationInterfaces(uSPCell,lagrangeP)				! -> matrice 4x(P+1)		(4 cotes et (P+1)

end extrapolationInterfaces


! uFace1 : inconnue sur la face 1
! uFace2 : inconnue sur la face 2
function uInterface(uFace1,uFace2)   						! -> vecteur (P+1)

end uInterface

! uCell : inconnue aux SP et aux interfaces d'une cellule (fonctionne aussi pour les flux?)
! lagrangeP : polynomes de Lagrange de degre P
! lagrangeP2 : polynomes de Lagrange de degre P+2
! dlagrangeP2 : derivees des polynomes de Lagrange de degre P+2
function naturalDerivatives(uCell,lagrangeP,lagrangeP2,dlagrangeP2) 		! -> matrice 2x(P+3)x(P+3)

end 

! Voir avis de G Puigt et structure des donnees de maillage ...
function realDerivatives()

end realDerivatives


! u : inconnue aux SP
! q : fonction auxiliaire (x) aux SP
! r : fonction auxiliaire (y) aux SP
! problemType = euler, navier-stokes,...
function flux(u,q,r,problemType)   						! -> matrice (P+1)x(P+1)

end uInterface

! uFace1 : inconnue sur la face 1
! uFace2 : inconnue sur la face 2
! fluxFace1 : flux sur la face 1
! fluxFace2 : flux sur la face 2
! tau : parametre de penalisation
function fluxInterface(uFace1,uFace2,fluxFace1,fluxFace2,tau)   						! -> vecteur 4x(P+1)

end uInterface

! u : inconnue aux SP
! dfdx : derivee du flux f selon x aux SP
! dgdy : derivee du flux g selon y aux SP
! alpha : parametres de RK
function timeMarch(u,dfdx,dgdy,alpha)				! -> vecteur (P+1)

end