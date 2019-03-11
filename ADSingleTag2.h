#ifndef ADSingleTag2_H
#define ADSingleTag2_H

//An example of making a tag using c++
//Based on UserTagging/UserTag/DetectorTag.py by Chao
//hem@ihep.ac.cn 2011-05-31
#include "Event/RecHeader.h"
#include "Event/JobInfo.h"
#include "Conventions/DetectorId.h"
#include "Conventions/JobId.h"
#include "Conventions/Trigger.h"
#include "DataSvc/IJobInfoSvc.h"
#include "DataUtilities/DybArchiveList.h"


#include "DybAlg/DybBaseAlg.h"
#include "Context/TimeStamp.h"
#include <vector>
#include <map>
#include <string>
#include <deque>

using namespace std;
class IJobInfoSvc;
class IDataProvider;
class TimeStamp;
namespace DayaBay
{
  class HeaderObject;
}
using namespace DayaBay;
typedef long long nstime_t;
typedef struct
{
    TimeStamp time;
    Detector detector;
    unsigned int triggerNumber;
} CacheItem;
class ADSingleTag2: public DybBaseAlg
{
  public:
    /// Constructor has to be in this form
    ADSingleTag2(const std::string& name, ISvcLocator* svcloc);
    virtual ~ADSingleTag2();

    /// Three mandatory member functions of any algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

  private:
    IJobInfoSvc *m_jobInfoSvc;
    IDataProviderSvc* p_archiveSvc;
    int timeCut;

    // Sam's methods
    bool isADMuon();
    bool isADShowerMuon();
    bool isWPMuon();
    bool isFlasher();
    bool hasPromptLikeEnergy();
    bool hasDelayedLikeEnergy();
    bool isVetoedByPastMuons();
    bool isVetoedByPastADSignal();
    bool isPromptLikeSingle();
    bool isDelayedLikeSingle();
    void maybeLogEventTimestamp();
    void initializeTimestampLogs();
    bool isVetoedByCurrent(TimeStamp cachedTime);
    bool isLaterThanAllVetoes(TimeStamp cachedTime);
    void checkCachesForVetoesAndPasses();
    StatusCode tagSpecifiedEvent(CacheItem event, string destination);
    static nstime_t timeStampToNS(const TimeStamp& t);
    // Sam's variables
    // convention: last_ refers to the past. cached_ refers to saved events
    // that might still pass the cuts. All other variables refer to the current
    // event.
    vector<unsigned int> numbersOfSingles;
    vector<unsigned int> numbersOfFailures;
    map<DetectorId::DetectorId_t, TimeStamp> last_AD_muon_time;
    map<DetectorId::DetectorId_t, TimeStamp> last_AD_shower_time;
    // The WP muon times apply to all detectors in a given EH (i.e. all in a
    // given file), so no map is needed.
    TimeStamp last_WP_muon_time;
    map<DetectorId::DetectorId_t, TimeStamp> last_IBD_like_time;
    deque<CacheItem> promptLikeCache;
    deque<CacheItem> delayedLikeCache;
    vector<deque<CacheItem>* > caches;
    vector<string> tagNames;
    bool didInitializeTimestampLogs;
    double Erec;
    double nominalCharge;
    int nHit;
    TimeStamp time;
    unsigned int triggerNumber;
    unsigned int triggerType;
    double MaxQRatio;
    double quadrantRatio;
    bool isVetoedByPast;
    Detector detector;
    static const nstime_t LAST_AD_MUON_VETO_NS = 1400000;
    static const nstime_t LAST_AD_SHOWER_VETO_NS = 400000000;
    static const nstime_t LAST_WP_MUON_VETO_NS = 600000;
    static const nstime_t LAST_IBD_LIKE_VETO_NS = 400000;
    static const nstime_t NEXT_AD_MUON_VETO_NS = 0;
    static const nstime_t NEXT_AD_SHOWER_VETO_NS = 0;
    static const nstime_t NEXT_WP_MUON_VETO_NS = 2000;
    static const nstime_t NEXT_IBD_LIKE_VETO_NS = 400000;
    static const double FLASHER_CONSTANT = 0.45;
    static const double FLASHER_THRESHOLD = 1.0;
    static const double PROMPT_ENERGY_MIN_MEV = 0.7;
    static const double PROMPT_ENERGY_MAX_MEV = 12.0;
    static const double DELAYED_ENERGY_MIN_MEV = 6.0;
    static const double DELAYED_ENERGY_MAX_MEV = 12.0;
    static const double AD_MUON_PE_THRESHOLD = 3e3;
    static const double AD_SHOWER_PE_THRESHOLD = 3e5;
    static const int WP_MUON_HIT_THRESHOLD = 12;
    //DetectorId::DetectorId_t did = recHeader->context().GetDetId();
    //TimeStamp currentTime = recHeader->timeStamp();
    //m_lastTimeStamp_map[did] = currentTime;


    //TimeStamp t;
    //DetectorId::DetectorId_t id;
    //t = m_lastTimeStamp_map[id];
};

#endif
