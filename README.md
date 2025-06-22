# Transcriptomics Casus
____
## 🧾 Inleiding 

Reumatoïde artritis (RA) is een langdurige auto-immuunziekte waarbij het afweersysteem het eigen lichaam aanvalt. Dit zorgt voor ontstekingen in het gewrichtsslijmvlies (synovium), wat uiteindelijk leidt tot schade aan de gewrichten (Smolen et al., 2016). Hoe RA precies ontstaat is nog niet helemaal duidelijk. Waarschijnlijk spelen erfelijke aanleg, omgevingsfactoren (zoals roken) en een verstoord immuunsysteem een rol (Firestein & McInnes, 2017). In het bloed van mensen met RA zijn vaak auto-antistoffen aanwezig, zoals anti-CCP (anti-cyclisch gecitrullineerd peptide). 
Deze helpen bij het stellen van een vroege diagnose (van der Helm-van Mil et al., 2025).

Het is belangrijk om RA vroeg te herkennen en te behandelen met medicijnen (DMARDs), zodat verdere gewrichtsschade wordt beperkt (Smolen et al., 2018). Tegenwoordig wordt met RNA-sequencing gekeken naar genexpressie in RA-weefsel. Dit geeft inzicht in welke genen en biologische processen bij RA betrokken zijn (Guo et al., 2023). In dit onderzoek worden RNA-gegevens geanalyseerd van synoviumbiopten van RA-patiënten en gezonde mensen. Het doel is om verschillen in genexpressie te vinden en te onderzoeken welke pathways mogelijk betrokken zijn bij het ontstaan van RA.

___
## 🧪 Methode

Voor dit onderzoek werd een RNA-sequencinganalyse uitgevoerd met monsters afkomstig van acht patiënten, vier hiervan met reumatoïde artritis (RA) (diagnoseduur > 12 maanden, anti-CCP-positief) en vier patiënten dienden als controlegroep zonder RA (anti-CCP-negatief). 

De ruwe data werd gedownload en ge-unzipt, vervolgens werd in RStudio gebruik gemaakt van de [`Rsubread`](./scripts/Rsubread.R) package om de referentie-index van het humane genoom (GRCh38.p14) (Homo Sapiens Genome Assembly GRCh38.p14, z.d.) op te bouwen en de reads te alignen naar het genoom. De uitgelijnde BAM-bestanden werden gesorteerd en geïndexeerd met [`Rsamtools`](./scripts/Rsamtools.R).

Met featureCounts werden de read counts per gen bepaald en hieruit volgde een count-matrix. Omdat dit nog een subset was, werd er vanaf nu gewerkt met de (aangereikte) volledige count-matrix. Er werd een differentiale genexpressieanalyse (DGE) uitgevoerd met de package [`DESeq2`](./scripts/deseq.R). 

Hierbij werd onder andere de log2-fold change (log2FC) en de adjusted p-value (padj) berekend, om de genen met significante expressieveranderingen tussen RA en de controlegroep eruit te halen. Resultaten werden weergegeven in o.a. een volcano plot (via [`EnhancedVolcanoplot`](./scripts/EnhancedVolcano.R)).

Vervolgens werd een Gene Ontology (GO) verrijkingsanalyse uitgevoerd met [`goseq`](./scripts/goseq.R) en visualisatie met ggplot2. Daarnast werd een KEGG-pathway onderzocht (op basis van de resultaten van de GO-analyse) met behulp van [`pathview`](./scripts/pathview.R), er werd gekozen voor pathway 'hsa04062', de pathway voor Chemokine signaling. Chemokines spelen een grote rol in het aantrekken van immuuncellen naar ontstekingsweefsel en bij RA is er chronische ontsteking van het synovium en daar zijn chemokines heel actief. Daarnaast reguleert deze pathway de migratie en activering van leukocyten (witte bloedcellen), deze leukocyt-gerelateerde processen kwamen ook naar voren in de [`GO-analyse`](./resultaten/GO-analyse.png) (400 counts, 32,8% hits en p-value = 3,847567e-22).

Scripts en data zijn te vinden in de mappen [`scripts`](./scripts) en [`data`](./data), het volledige script is [`Rscript`](./scripts/Rscript.R)

<p align="center">
  <img src="assets/Workflow.png" alt="Workflow RNA-Seq data-analyse bij RA" width="600"/>
</p>

___
## 📊 Resultaten

Om de verschillen in genexpressie te vinden tussen patiënten met reumatoïde artritis (RA) en gezonde controles, is een differentiële genexpressie-analyse uitgevoerd. In de [`volcano plot`](./resultaten/EnhancedVolcanoPlot.png) zijn de genen te zien die significant verschillend tot expressie komen. De genen in het rood hebben zowel een statistische significante p-waarde als een hoge log₂ fold change, met deze genen wordt dan ook verder gewerkt in de analyse.

Hierna werd een [`GO (Gene Ontology) analyse`](./resultaten/GO-analyse.png) uitgevoerd en hieruit kwam dat immuungerelateerde processen meer betrokken zijn bij patiënten met RA. Dit sluit aan bij RA als een auto-immuunziekte, waarbij dus de immuunactivatie een grote rol speelt. Zo is ook te zien dat "leukocyte activation" in de GO-analyse naar voren komt  met 400 counts, 32,8% hits en een p-value = 3,847567e-22. Met dit proces is verder gezocht naar pathways voor de KEGG-pathway analyse.

Voor deze KEGG-pathway analyse werd de [`chemokine signaling pathway`](./resultaten/hsa04062.pathview.png) geanalyseerd. Binnen deze pathway is te zien dat er meerdere genen differentieel gereguleerd zijn tussen RA-patiënten en gezonde individuen. Onder andere PI3K, PLC, Rac en PKC betrokken bij processen als leukocytenmigratie, ontsteking en celactivatie. De rode kleur staat voor upregulatie, wat betekent dat het gen meer tot expressie is gekomen in de RA-patiënten dan in gezonde patiënten (Een Log₂FC = +5 betekent ongeveer 32x meer tot expressie). Groen staat juist voor downregulatie, waarbij het gen minder tot expressie is gekomen in RA-patiënten dan in de controles (Een Log₂FC = -5 betekent ongeveer 32x minder tot expressie). De activatie van deze route laat de rol zien van chemokines in het immuunantwoord bij RA. 
___
## 🎯 Conclusie


___
## 🗂️ Data Stewardship

- Zie [`data_stewardship`](./data_stewardship) voor uitleg over het beheren van data volgens de FAIR-principes.

___

*Gemaakt door Rika Ferwerda*
