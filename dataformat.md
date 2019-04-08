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
 - NominalCharge (number of photoelectrons)
 - Maximum charge on any 2-inch PMT
 - Flasher PSD t1 (200ns/400ns)
 - Flasher PSD t2 (150ns/400ns)
 - fQuad = Q3/(Q2+Q4)
 - fMax = Q(highest charge)/Qtotal

## Quantities to compute and then save

### Tags

 - Flasher
 - 2-inch flasher
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

 - [X] flasher ID fID
 - [X] flasher PSD ID fPSD
 - [X] number of muons within previous 5s
 - [X] {time to each muon within previous 5s}
 - [X] time to previous WP muon
 - [X] time to next WP muon
 - [X] time to previous AD muon
 - [] time to next AD muon
 - [X] time to previous AD shower muon
 - [] time to next AD shower muon
 - [X] number of prompt-like signals in previous 400us
 - [X] {time to each prompt-like signal in previous 400us}
 - [X] time to previous prompt-like signal
 - [] time to next delayed-like signal
