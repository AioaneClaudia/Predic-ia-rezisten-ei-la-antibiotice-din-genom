# Raport pilot: date BV-BRC pentru E. coli

## Date folosite
- Fenotipurile de rezistenta si metadatele genomurilor au fost descarcate din BV-BRC prin API (`src/01_descarca_date_api.py`).
- Am ales doua antibiotice: ciprofloxacin (rezistenta prin mutatii in cromozom) si ampicilina (gene mobile).
- Din 489.375 de randuri AMR, aproximativ 460.000 sunt predictii ale BV-BRC (modele XGBoost). Le-am exclus si am pastrat doar rezultatele de laborator (29.845 de randuri).

## Curatarea datelor
- Acelasi genom apare uneori de mai multe ori. Sunt 15.675 de perechi genom x antibiotic, dintre care 1.523 cu mai multe inregistrari.
- Doar 41 de perechi au rezultate contradictorii. Le-am eliminat, la fel ca valorile intermediare (154).
- Au ramas 8.561 de perechi sensibile si 6.919 rezistente.
- Genome unice: 6.364 pentru ciprofloxacin si 5.507 pentru ampicilina.

## Anul de izolare
- Anul exista pentru majoritatea genomelor. Intre 2001 si 2017 avem sute de izolate pe an.
- Dupa 2017 datele sunt putine si instabile (rezistenta la ciprofloxacin sare de la 13% in 2021 la 85% in 2022), deci ne oprim la 2017.

## Sursa datelor
- Norvegia are cele mai multe genome (3.156 pentru ciprofloxacin), din 2002 pana in 2017. Numarul si perioada seamana cu un studiu national de supraveghere a infectiilor sanguine cu E. coli. Nu e confirmat, pentru ca aproape toate genomele norvegiene nu au bioproiect completat.
- Marea Britanie are 1.489 de genome pentru ciprofloxacin, aproape toate din proiectul PRJEB4681.
- In datele combinate, rezistenta la ciprofloxacin pare sa creasca brusc dupa 2012. In seria norvegiana singura, cresterea este treptata (de la 4% in 2002 la 15% in 2017). Deci saltul vine din amestecul surselor.

## Decizie pentru proof of concept
- Folosim ciprofloxacin si doar genomele din Norvegia, 2002-2017, ca sa nu confundam timpul cu tara sau cu studiul.
- Impartire temporala: antrenare 2002-2012 (114 rezistente, 1.612 sensibile) si test 2013-2017 (195 rezistente, 1.230 sensibile).
- Prevalenta rezistentei creste de la aproximativ 7% in antrenare la aproximativ 14% in test, deci folosim AUPRC si MCC, nu acuratetea.
- Toate genomele norvegiene au ST si clustere cgMLST (HC0-HC100) in metadate, iar 99,9% au calitate buna. Le folosim pentru clustere filogenetice si pentru baseline-ul doar-linie, fara sa descarcam genome.

## Riscuri
- Sursa norvegiana nu este confirmata de bioproiect.
- Avem doar 309 izolate rezistente in total (aproximativ 10%).
- Semnificatia coloanelor `cgmlst_hc*` trebuie verificata in documentatia BV-BRC.

## Urmatorii pasi
- Notebook 02: baseline doar-linie cu trei scheme de impartire (aleatorie, temporala, pe clustere).
- Masurarea distantei dintre test si antrenare cu nivelurile HC.
- Descarcarea genomelor si modelul pe k-meri.