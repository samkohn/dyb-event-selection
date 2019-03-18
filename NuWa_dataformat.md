My description of the NuWa output data format
=============================================

## recon.Neutrino.\*.root files

Contain calibrated + reconstructed event info.

### Event/Data/CalibStats

Metadata:
- Event timestamp (context.mTimeStamp)
- EH (context.mSite)
- Detector ID (context.mDetId)
- Trigger number (triggerNumber)

Physics data:
- Trigger time (context.mTimeStamp)
- NHit (nHit)
- fQuad (Quadrant)
- fMax (MaxQ)
- Flasher PSD t1 (time\_PSD)
- Flasher PSD t2 (time\_PSD1)
- Maximum charge for 2-inch PMTs (MaxQ\_2inchPMT)

### Event/Rec/AdSimple

Metadata:
- Trigger type (triggerType)

Physics data:
- Energy (energy)
- Reconstructed position (x, y, z)

