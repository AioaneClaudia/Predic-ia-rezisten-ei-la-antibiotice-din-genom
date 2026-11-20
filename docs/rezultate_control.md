# Rezultate: controlul cu distanta potrivita (ciprofloxacin, Norvegia 2002-2017)

## Setup
- Model doar-linie la nivel HC50. Antrenare: 65% din izolatele din 2002-2012 (circa 1.119); test temporal: 2013-2017 (1.425 izolate).
- Control: izolate din 2002-2012 neutilizate la antrenare, esantionate (circa 400) astfel incat distributia distantei pana la antrenare sa fie aceeasi ca la testul temporal (verificat: de exemplu HC50 40,4% in ambele).
- 300 de repetari, cu bootstrap pe izolatele de test.

## Rezultate
1. Impartiri pe clustere: modelele pe un nivel HC cad la AUROC 0,5 cand clusterele acelui nivel sunt separate intre antrenare si test; modelele pe niveluri mai largi isi pastreaza informatia (HC100 la split HC50: 0,732).
2. Control: AUROC temporal 0,843 [0,805; 0,878], AUROC control 0,867 [0,788; 0,940]; diferenta -0,024 [-0,104; 0,065]. Intervalul include 0.
3. Controlul are circa 25 de izolate rezistente, deci intervalele sunt largi; castigul AUPRC nu poate fi interpretat (diferenta 0,027, interval [-0,227; 0,301]).
4. Distributia distantei pana la antrenare nu se schimba in anii 2013-2017 (de exemplu HC50: 40-43% in fiecare an), deci distanta apare ca un salt la momentul taierii.
5. Regresie la nivel de izolat: efectul decalajului asupra erorii este -0,0021 pe an [-0,0091; 0,0049], p = 0,56, si ramane -0,0022 dupa includerea distantei; pentru reziduali -0,0013 si -0,0015 (p = 0,73-0,76).

## Interpretare prudenta
Pentru modelul doar-linie nu detectam un efect al timpului peste cel al distantei si nu vedem scadere progresiva in anii 2013-2017. Puterea statistica este limitata (circa 25 de rezistente in control), deci rezultatul inseamna absenta unui efect detectabil, nu absenta unui efect.

## Limite
- Valabil doar pentru modelul doar-linie; urmeaza modelele pe k-meri.
- Impartirea pe clustere la HC100 a produs un test de 55% din date (clustere mari).
- Un singur antibiotic si o singura sursa; ampicilina si seria din Marea Britanie urmeaza.

## Comparatie ciprofloxacin si ampicilina (model doar-linie)
- Linia explica mult mai bine rezistenta la ciprofloxacin (AUROC temporal 0,877 cu ST, 0,897 cu HC100) decat la ampicilina (0,689 si 0,726). Este in acord cu ideea ca rezistenta mobila depinde mai putin de linie.
- Ciprofloxacin: diferenta temporal - rezervor fara potrivire = -0,053 [-0,117; 0,019]; cu potrivire pe distanta = -0,024 [-0,104; 0,065]. Aproximativ jumatate din scaderea aparenta se explica prin distanta, dar intervalele includ 0.
- Ampicilina: diferenta fara potrivire = 0,003 [-0,059; 0,067]; cu potrivire = 0,031 [-0,042; 0,108]. Nu exista o scadere de explicat.
- Regresia la nivel de izolat: niciun efect semnificativ al decalajului temporal (ciprofloxacin: -0,0021 pe an, p = 0,56; ampicilina: +0,0033 pe an, p = 0,33), nemodificat de includerea distantei.
- Controlul are circa 25 de rezistente la ciprofloxacin si circa 138 la ampicilina; intervalele raman largi.

## Limite
- Valabil doar pentru modelul doar-linie; tipul de distanta (nucleu fata de accesoriu) nu a fost testat.