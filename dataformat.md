Data format notes
=================

## Quantities to save from basic reconstructed data

### Metadata

 - Run number
 - Run start date/time
 - Run end date/time
 - EH
 - Detector (AD/WP/etc)
 - Trigger number
 - Trigger type

### Physics data

 - Energy
 - Reconstructed position
 - Trigger time
 - NHit (number of hit PMTs)
 - Maximum charge on any 2-inch PMT
 - Flasher PSD t1 (200ns/400ns)
 - Flasher PSD t2 (150ns/400ns)
 - fQuad = Q3/(Q2+Q4)
 - fMax = Q(highest charge)/Qtotal

## Quantities to compute and then save

### Tags

 - Muon
   - WP Muon
   - AD Muon
   - AD Shower Muon
 - Single events
   - Prompt-like
   - Delayed-like
 - Fast neutron events
   - Boundary-muon-tagged low-energy prompt
   - Boundary-muon-tagged low-energy delayed
   - Boundary-muon-tagged high-energy prompt
   - Boundary-muon-tagged high-energy delayed
   - IBD-like high energy prompt
   - IBD-like high energy delayed
 - Li9/He8/B12
 - Am-C
   - intense-run single-event
   - intense-run correlated prompt
   - intense-run correlated delayed
   - regular single-event

### Processed data

 - number of muons within previous 5s
 - [time to each muon within previous 5s]
 - time to previous WP muon
 - time to next WP muon
 - time to previous AD muon
 - time to next AD muon
 - time to previous AD shower muon
 - time to next AD shower muon
