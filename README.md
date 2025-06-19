# Transcriptomics Casus
____
## 🧾 Inleiding 

Reumatoïde artritis (RA) is een langdurige auto-immuunziekte waarbij het afweersysteem het eigen lichaam aanvalt. Dit zorgt voor ontstekingen in het gewrichtsslijmvlies (synovium), wat uiteindelijk leidt tot schade aan de gewrichten (Smolen et al., 2016). Hoe RA precies ontstaat is nog niet helemaal duidelijk. Waarschijnlijk spelen erfelijke aanleg, omgevingsfactoren (zoals roken) en een verstoord immuunsysteem een rol (Firestein & McInnes, 2017). In het bloed van mensen met RA zijn vaak auto-antistoffen aanwezig, zoals anti-CCP (anti-cyclisch gecitrullineerd peptide). 
Deze helpen bij het stellen van een vroege diagnose (van der Helm-van Mil et al., 2025).

Het is belangrijk om RA vroeg te herkennen en te behandelen met medicijnen (DMARDs), zodat verdere gewrichtsschade wordt beperkt (Smolen et al., 2018). Tegenwoordig wordt met RNA-sequencing gekeken naar genexpressie in RA-weefsel. Dit geeft inzicht in welke genen en biologische processen bij RA betrokken zijn (Guo et al., 2023). In dit project worden RNA-gegevens geanalyseerd van synoviumbiopten van RA-patiënten en gezonde mensen. Het doel is om verschillen in genexpressie te vinden en te onderzoeken welke pathways mogelijk betrokken zijn bij het ontstaan van RA.

___
## 🧪 Methode

Voor dit onderzoek werd een RNA-sequencinganalyse uitgevoerd met monsters afkomstig van acht patiënten, vier hiervan met reumatoïde artritis (RA) (diagnoseduur > 12 maanden, anti-CCP-positief) en vier patiënten dienden als controlegroep zonder RA (anti-CCP-negatief). 

De ruwe data werd gedownload en ge-unzipt, vervolgens werd in RStudio gebruik gemaakt van de Rsubread-package om de referentie-index van het humane genoom (GRCh38.p14) op te bouwen en de reads te alignen naar het genoom. De uitgelijnde BAM-bestanden werden gesorteerd en geïndexeerd met Rsamtools.
Met featureCounts werden de read counts per gen bepaald en hieruit volgde een count-matrix. Omdat dit nog een subset was, werd er vanaf nu gewerkt met de (aangereikte) volledige count-matrix. Er werd een differentiale genexpressieanalyse (DGE) uitgevoerd met het package DESeq2. Hierbij werd onder andere de log2-fold change (log2FC) en de adjusted p-value (padj) berekend, om de genen met significante expressieveranderingen tussen RA en de controlegroep eruit te halen.
Vervolgens werd een Gene Ontology (GO) verrijkingsanalyse uitgevoerd met goseq en visualisatie met ggplot2. Daarnast werden KEGG-pathways onderzocht met behulp van pathview. Resultaten werden weergegeven in o.a. een volcano plot (via EnhancedVolcano).

Scripts en data zijn te vinden in de mappen [`scripts`](./scripts) en [`data`](./data)  
___
## 📊 Resultaten


___
## 🎯 Conclusie


___
## 🗂️ Data Stewardship

- Zie [`data_stewardship`](./data_stewardship) voor uitleg over het beheren van data volgens de FAIR-principes.

___

*Gemaakt door Rika Ferwerda*
