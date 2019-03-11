#include "ADSingleTag2.h"
#include "Event/ReadoutHeader.h"
#include "Event/CalibReadoutHeader.h"
#include "Event/UserDataHeader.h"
#include "Event/JobInfo.h"
#include "Conventions/DetectorId.h"
#include "Conventions/JobId.h"
#include "DataSvc/IJobInfoSvc.h"
#include "Context/ServiceMode.h"
#include "Event/RecHeader.h"
#include "Event/RecTrigger.h"
#include <math.h>
#include "DataUtilities/DybArchiveList.h"

using namespace std;
using namespace DayaBay;

ADSingleTag2::ADSingleTag2(const string& name, ISvcLocator* svcloc) :
    DybBaseAlg(name, svcloc) {
        // this is just an example
        declareProperty("dtCut",timeCut,
                "single time Cut");
    }

ADSingleTag2::~ADSingleTag2() {
}

StatusCode ADSingleTag2::initialize() {
    debug() << "initialize()" << endreq;

    m_jobInfoSvc = svc<IJobInfoSvc> ("JobInfoSvc", true);
    if (!m_jobInfoSvc) {
        error() << "Failed to initialize JobInfoSvc" << endreq;
        return StatusCode::FAILURE;
    }
    //Get Archive Svc
    StatusCode status = service("EventDataArchiveSvc", p_archiveSvc);
    if (status.isFailure()) {
        Error("Service [EventDataArchiveSvc] not found", status);
    }

    // I would like to initialize the last_<...>_time veto records to the
    // timestamp of the first event here but I do not have access to that until
    // the execute() method, so I will initialize there.
    didInitializeTimestampLogs = false;

    // Set up the caches
    caches.push_back(&promptLikeCache);
    caches.push_back(&delayedLikeCache);

    // Set up the tag names
    tagNames.push_back("/Event/Tag/Kohn/SinglesPromptLike");
    tagNames.push_back("/Event/Tag/Kohn/SinglesDelayedLike");

    // Set up the singles counters (with the same number order as the
    // caches and tagnames vectors).
    numbersOfSingles.push_back(0);
    numbersOfSingles.push_back(0);
    numbersOfFailures.push_back(0);
    numbersOfFailures.push_back(0);

    return StatusCode::SUCCESS;
}

StatusCode ADSingleTag2::execute() {
    /*
     * This method has 5 sections:
     *   - Fetch data from other data structures
     *   - Compute important quantities and calculate vetoes
     *   - Update the log of past muons/other signals
     *   - Determine if any cached events are vetoed or accepted
     *   - Decide whether to cache this event
     */
    debug() << "execute() ______________________________ start" << endreq;

    // ****************************************************

    // Fetch all of the data for the current event
    // Some values are derived in reconstruction. Some are in CalibStats.
    // Here I get the NuWa objects (RecHeader -> reconstruction and
    // calibStatsHeader -> CalibStats) required for the data.
    RecHeader* recHeader = get<RecHeader>("/Event/Rec/AdSimple");
    const RecTrigger& recTrigger = recHeader->recTrigger();
    UserDataHeader* calibStatsHeader = 0;
    calibStatsHeader = get<UserDataHeader>("/Event/Data/CalibStats");

    // Store the data in this object.
    this->time = recTrigger.triggerTime();
    this->triggerNumber = recTrigger.triggerNumber();
    this->Erec = recTrigger.energy();
    this->detector = recTrigger.detector();
    this->nominalCharge = calibStatsHeader->getFloat("NominalCharge");
    this->nHit = calibStatsHeader->getInt("nHit");
    this->MaxQRatio = calibStatsHeader->getFloat("MaxQ");
    this->quadrantRatio = calibStatsHeader->getFloat("Quadrant");

    // ****************************************************

    // If the current event is a flasher, immediately skip it without updating
    // any logs or doing anything else. Also tag it as a flasher for further
    // analysis
    if(isFlasher())
    {
        // Only save 1/8 of the flashers
        if(time.GetNanoSec() % 8 == 0)
        {
            string destination = "/Event/Tag/Kohn/Flasher";
            CalibReadoutHeader* readoutHeader = get<CalibReadoutHeader>("/Event/CalibReadout/CalibReadoutHeader");
            vector<const IHeader*> inputHeaders;
            inputHeaders.push_back(calibStatsHeader);
            inputHeaders.push_back(readoutHeader);
            inputHeaders.push_back(recHeader);
            HeaderObject* singlesTag = MakeHeader<HeaderObject>(inputHeaders);
            debug() << "Saving flasher trigger number " << readoutHeader->calibReadout()->triggerNumber() << endreq;
            put(singlesTag, destination.c_str());
        }
        return StatusCode::SUCCESS;
    }

    // Determine multiplicity vetoes
    // 0th task is to initialize timestamp logs if not already done
    if(!didInitializeTimestampLogs)
    {
        initializeTimestampLogs();
    }

    // Only check for vetoes if event happened in an AD (otherwise return/exit
    // in a few lines). This check is performed before the current event is
    // logged so that the event will not veto itself.
    if(detector.isAD())
    {
        isVetoedByPast = isVetoedByPastMuons() || isVetoedByPastADSignal();
    }

    // ****************************************************
    maybeLogEventTimestamp();

    // ****************************************************
    checkCachesForVetoesAndPasses();

    // ****************************************************

    // Here we return without tagging if the event happened not in an AD.
    //if(!detector.isAD())
    //{
        //return StatusCode::SUCCESS;
    //}

    if(isPromptLikeSingle())
    {
        CacheItem event;
        event.time = time;
        event.detector = detector;
        event.triggerNumber = triggerNumber;
        promptLikeCache.push_back(event);
    }
    if(isDelayedLikeSingle())
    {
        CacheItem event;
        event.time = time;
        event.detector = detector;
        event.triggerNumber = triggerNumber;
        delayedLikeCache.push_back(event);
    }

debug() << "execute() ______________________________ end" << endreq;
return StatusCode::SUCCESS;

}

StatusCode ADSingleTag2::finalize() {
    //  string deadtimefile="test.txt";
    //    ofstream outputfile;
    //      outputfile.open(deadtimefile.c_str());
    debug() << "ADSingleTag2::finalize()" << endreq;
    unsigned int numberOfPrompts = numbersOfSingles.at(0);
    unsigned int numberOfDelayeds = numbersOfSingles.at(1);
    info() << "ADSingleTag2 found " << numberOfPrompts << " prompt-like singles and ";
    info() << numberOfDelayeds << " delayed-like singles." << endreq;
    info() << "ADSingleTag2 failed " << numbersOfFailures.at(0) + numbersOfFailures.at(1);
    info() << " times." << endreq;
    //  outputfile<< deadtime<< endl;
    return StatusCode::SUCCESS;
}

bool ADSingleTag2::hasPromptLikeEnergy()
{
    return (Erec > PROMPT_ENERGY_MIN_MEV &&
            Erec < PROMPT_ENERGY_MAX_MEV);
}

bool ADSingleTag2::hasDelayedLikeEnergy()
{
    return (Erec > DELAYED_ENERGY_MIN_MEV &&
            Erec < DELAYED_ENERGY_MAX_MEV);
}

/*
 * The criterion to pass is (not a muon) and (has prompt-like energy)
 * and (not vetoed by past muons) and (not vetoed by past AD triggers)
 * and (detector is an AD).
 */
bool ADSingleTag2::isPromptLikeSingle()
{
    return (!isADMuon() &&
            !isADShowerMuon() &&
            !isWPMuon() &&
            hasPromptLikeEnergy() &&
            !isVetoedByPast &&
            detector.isAD());
}

/*
 * The criterion to pass is (not a muon) and (has delayed-like energy)
 * and (not vetoed by past muons) and (not vetoed by past AD triggers)
 * and (detector is an AD).
 */
bool ADSingleTag2::isDelayedLikeSingle()
{
    return (!isADMuon() &&
            !isADShowerMuon() &&
            !isWPMuon() &&
            hasDelayedLikeEnergy() &&
            !isVetoedByPast &&
            detector.isAD());
}

/*
 * AD Muons have >= 3000 photoelectrons (nominal charge) and are located in an
 * AD (detectors 1 through 4).
 */
bool ADSingleTag2::isADMuon()
{
    bool isAD = detector.isAD();
    bool hasEnoughPhotoelectrons = (nominalCharge >= AD_MUON_PE_THRESHOLD);
    return (hasEnoughPhotoelectrons && isAD);
}

/*
 * AD Shower Muons have >= 30000 photoelectrons (nominal charge) and are
 * located in an AD.
 */
bool ADSingleTag2::isADShowerMuon()
{
    bool isAD = detector.isAD();
    bool hasEnoughPhotoelectrons = (nominalCharge >= AD_SHOWER_PE_THRESHOLD);
    return (hasEnoughPhotoelectrons && isAD);
}

/*
 * Water Pool Muons have NHIT > 12 (>12 PMTs hit) and are located in a water
 * pool.
 */
bool ADSingleTag2::isWPMuon()
{
    bool isWP = detector.isWaterShield();
    bool hasEnoughHits = nHit > WP_MUON_HIT_THRESHOLD;
    return (hasEnoughHits && isWP);
}

/*
 * Flashers follow the formula specified in the code below.
 */
bool ADSingleTag2::isFlasher()
{
    double scaledQuadrant = quadrantRatio/FLASHER_CONSTANT;
    double flasherParam = MaxQRatio * MaxQRatio + scaledQuadrant * scaledQuadrant;
    return flasherParam > FLASHER_THRESHOLD;
}

/*
 * Conditionally log the event timestamp if it has the potential to veto
 * future events (e.g. if it is a muon).
 */
void ADSingleTag2::maybeLogEventTimestamp()
{
    bool isMuon = false;
    if(isADMuon())
    {
        last_AD_muon_time[detector.detectorId()] = time;
        isMuon = true;
    }
    if(isADShowerMuon())
    {
        last_AD_shower_time[detector.detectorId()] = time;
        isMuon = true;
    }
    if(isWPMuon())
    {
        last_WP_muon_time = time;
        isMuon = true;
    }
    if(!isMuon && (hasPromptLikeEnergy() || hasDelayedLikeEnergy()))
    {
        last_IBD_like_time[detector.detectorId()] = time;
    }
}

/*
 * Compare the last muon times to the thresholds given as static consts
 * in this class. A return value of true indicates the event should be vetoed.
 *
 * The initial conditions for the first few events assume that event 0
 * is the most recent event of all veto-able types.
 */
bool ADSingleTag2::isVetoedByPastMuons()
{
    TimeStamp dt_AD_muon = time - last_AD_muon_time[detector.detectorId()];
    TimeStamp dt_AD_shower = time - last_AD_shower_time[detector.detectorId()];
    TimeStamp dt_WP_muon = time - last_WP_muon_time;
    // nstime_t can go up to +/- 9e18, which should be enough for the
    // next 100 or so years (we're currently at 1.4e18 ns since the epoch).
    nstime_t dt_AD_muon_ns = timeStampToNS(dt_AD_muon);
    nstime_t dt_AD_shower_ns = timeStampToNS(dt_AD_shower);
    nstime_t dt_WP_muon_ns = timeStampToNS(dt_WP_muon);

    // Check each dt for being positive and veto if not (since that
    // means the events are out of order and the time veto is not
    // reliable).
    if(dt_AD_muon_ns < 0 ||
            dt_AD_shower_ns < 0 ||
            dt_WP_muon_ns < 0)
    {
        return true;
    }
    // Check each dt against the relevant threshold
    bool isVetoedByADMuon = dt_AD_muon_ns < LAST_AD_MUON_VETO_NS;
    bool isVetoedByADShower = dt_AD_shower_ns < LAST_AD_SHOWER_VETO_NS;
    bool isVetoedByWPMuon = dt_WP_muon_ns < LAST_WP_MUON_VETO_NS;

    return isVetoedByADMuon || isVetoedByADShower || isVetoedByWPMuon;
}

/*
 * Compare the last AD IBD-like times to the thresholds given as static
 * consts in this class. A reteurn value of true indicates the event
 * should be vetoed
 *
 * The initial conditions for the first few events assume that event 0
 * is the most recent event of all veto-able types.
 */
bool ADSingleTag2::isVetoedByPastADSignal()
{
    TimeStamp dt_AD_signal = time - last_IBD_like_time[detector.detectorId()];
    nstime_t dt_AD_signal_ns = timeStampToNS(dt_AD_signal);

    // Check each dt for being positive and veto if not (since that
    // means the events are out of order and the time veto is not
    // reliable).
    if(dt_AD_signal_ns < 0)
    {
        return true;
    }

    bool isVetoedByIBDLikeSignal = dt_AD_signal_ns < LAST_IBD_LIKE_VETO_NS;
    return isVetoedByIBDLikeSignal;
}

/*
 * Conservatively assume that the 0th event is "bad"/veto-worthy in all
 * possible ways, so that a long gap is required until the first accepted
 * event.
 */
void ADSingleTag2::initializeTimestampLogs()
{
    // When looping through detectors, kAll is the last one defined in the
    // enum, so it gets to be the "end."
    DetectorId::DetectorId_t end = DetectorId::kAll;
    // Loop through all DetectorId_t's using the duality between enums and
    // ints.
    for(int d = 0; d <= end; ++d)
    {
        DetectorId::DetectorId_t det = (DetectorId::DetectorId_t) d;
        last_AD_muon_time[det] = time;
        last_AD_shower_time[det] = time;
        last_IBD_like_time[det] = time;
    }
    last_WP_muon_time = time;
    didInitializeTimestampLogs = true;
    return;
}

/*
 * Cycle through the event caches and check for multiplicity vetoes.
 * Additionally, tag (finally) any cached events which survived the veto
 * process.
 *
 * The event caches are double-ended queues (deques) which have
 * efficient adding and deleting from both the front and back ends. The
 * back is the oldest event, which might finally be old enough that it
 * cannot be vetoed by any new event (so then we tag it and remove it
 * from the deque). The front is the most recent non-current event,
 * which is the most at-risk for being vetoed by the current event (in
 * which case we drop it).
 */
void ADSingleTag2::checkCachesForVetoesAndPasses()
{
    // Check all caches (e.g. prompt, delayed)
    for(size_t i = 0; i < caches.size(); ++i)
    {
        deque<CacheItem>& cache = *(caches.at(i));
        string tagName = tagNames.at(i);
        // Keep track of whether the cache changes. If it does not, then
        // keep_going should be false and the loop can stop.
        bool keep_going = true;
        while(keep_going && cache.size() > 0)
        {
            keep_going = false;
            // If cache.size() == 1 then oldest == newest. This would be
            // a problem if it got vetoedByCurrent and then was also
            // laterThanAllVetoes, but that is not logically possible.
            CacheItem oldest = cache.front();
            CacheItem newest = cache.back(); // Newest except for current event
            if(isVetoedByCurrent(newest.time))
            {
                debug() << "Cached event vetoed" << endreq;
                cache.pop_back();
                keep_going = true;
            }
            if(isLaterThanAllVetoes(oldest.time))
            {
                debug() << "Cached event accepted" << endreq;
                StatusCode status = tagSpecifiedEvent(oldest, tagName);
                cache.pop_front();
                if(status == StatusCode::SUCCESS)
                {
                    ++numbersOfSingles.at(i);
                }
                else
                {
                    ++numbersOfFailures.at(i);
                }
                keep_going = true;
            }
        }
    }
    return;
}

/*
 * Go through the Archive Event Service (AES) to find the specified cached
 * event and put() it to the given destination. Sometimes there are 2 TTree
 * events with the same TimeStamp (e.g. WP and AD). I believe there are never 2
 * AD events which share a TimeStamp. Hence, I will enforce a further cut that
 * the retrieved event must be from an AD.
 */
StatusCode ADSingleTag2::tagSpecifiedEvent(CacheItem event, string destination)
{
    SmartDataPtr<DybArchiveList> calibStatsList(p_archiveSvc, "/Event/Data/CalibStats");
    SmartDataPtr<DybArchiveList> calibReadoutList(p_archiveSvc, "/Event/CalibReadout/CalibReadoutHeader");
    SmartDataPtr<DybArchiveList> adSimpleList(p_archiveSvc, "/Event/Rec/AdSimple");
    if(!calibStatsList || !calibReadoutList || !adSimpleList)
    {
        error() << "Unable to fetch required headers from AES!!" << endreq;
        return StatusCode::FAILURE;
    }

    DybArchiveList::const_iterator statsIter = calibStatsList->begin();
    DybArchiveList::const_iterator readoutIter = calibReadoutList->begin();
    DybArchiveList::const_iterator recIter = adSimpleList->begin();

    UserDataHeader* statsHeaderToSave = 0;
    CalibReadoutHeader* readoutHeaderToSave = 0;
    RecHeader* recHeaderToSave = 0;
    bool keep_going = true;
    for(/*initialized above*/; keep_going && statsIter != calibStatsList->end(); ++statsIter, ++readoutIter, ++recIter)
    {
        UserDataHeader* statsHeader = dynamic_cast<UserDataHeader*> (*statsIter);
        CalibReadoutHeader* readoutHeader = dynamic_cast<CalibReadoutHeader*> (*readoutIter);
        RecHeader* recHeader = dynamic_cast<RecHeader*> (*recIter);

        // Check the timestamp to see if they're equal
        const CalibReadout* readout = readoutHeader->calibReadout();
        bool triggerTimeMatches = readout->triggerTime() == event.time;
        bool detectorMatches = readout->detector() == event.detector;
        bool triggerNumberMatches = readout->triggerNumber() == event.triggerNumber;
        if(triggerTimeMatches && triggerNumberMatches && detectorMatches)
        {
            // save and quit
            statsHeaderToSave = statsHeader;
            readoutHeaderToSave = readoutHeader;
            recHeaderToSave = recHeader;
            keep_going = false;
        }
    }
    if(keep_going)
    {
        return StatusCode::FAILURE;
    }
    vector<const IHeader*> inputHeaders;
    inputHeaders.push_back(statsHeaderToSave);
    inputHeaders.push_back(readoutHeaderToSave);
    inputHeaders.push_back(recHeaderToSave);
    HeaderObject* singlesTag = MakeHeader<HeaderObject>(inputHeaders);
    debug() << "Saving cached trigger number " << readoutHeaderToSave->calibReadout()->triggerNumber() << endreq;
    put(singlesTag, destination.c_str());
    return StatusCode::SUCCESS;
}

/*
 * Check to see if the given cached event is vetoed (mulitplicity-wise)
 * by the current event. Returns true if the cached event should be
 * discarded.
 */
bool ADSingleTag2::isVetoedByCurrent(TimeStamp cachedTime)
{
    TimeStamp dt = time - cachedTime;
    nstime_t dt_ns = timeStampToNS(dt);
    bool veto = false;
    if(isADMuon())
    {
        veto = veto || dt_ns < NEXT_AD_MUON_VETO_NS;
    }
    if(isADShowerMuon())
    {
        veto = veto || dt_ns < NEXT_AD_SHOWER_VETO_NS;
    }
    if(isWPMuon())
    {
        veto = veto || dt_ns < NEXT_WP_MUON_VETO_NS;
    }
    if(isPromptLikeSingle())
    {
        veto = veto || dt_ns < NEXT_IBD_LIKE_VETO_NS;
    }
    if(isDelayedLikeSingle())
    {
        veto = veto || dt_ns < NEXT_IBD_LIKE_VETO_NS;
    }
    return veto;
}

/*
 * Check to see if the given cached event is far enough in the past that
 * it is not possible to veto it anymore. This code might become invalid if the
 * veto constants change and there is a new longest veto time.
 */
bool ADSingleTag2::isLaterThanAllVetoes(TimeStamp cachedTime)
{
    nstime_t longest_wait = NEXT_IBD_LIKE_VETO_NS;
    TimeStamp dt = time - cachedTime;
    nstime_t dt_ns = timeStampToNS(dt);

    return dt_ns > longest_wait;
}

/*static*/ nstime_t ADSingleTag2::timeStampToNS(const TimeStamp& t)
{
    nstime_t s_to_ns = 1000000000;
    return s_to_ns * t.GetSec() + t.GetNanoSec();
}
