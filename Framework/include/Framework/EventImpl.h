/**
 * @file EventImpl.h
 * @brief Class implementing an event buffer system for storing event data
 * @author Jeremy McCormick, SLAC National Accelerator Laboratory
 */

#ifndef EVENT_EVENTIMPL_H_
#define EVENT_EVENTIMPL_H_

// ROOT
#include "TObject.h"
#include "TClonesArray.h"

// LDMX-SW
#include "Event/Event.h"

class TTree;
class TBranch;

// STL
#include <string>
#include <map>

namespace ldmxsw {

/**
 * @class EventImpl
 * @brief Implements an event buffer system for storing event data
 *
 * @note
 * Event data is stored in ROOT trees and branches, which can be added
 * on the fly.  The same TClonesArray and TObject pointers should be
 * used to add objects and collections from user code, as the class will
 * add a data structure or new ones automatically.
 */
  class EventImpl : public event::Event {

    public:

        /**
         * Class constructor.
         * @param passName The default pass name for adding event data.
         */
        EventImpl(const std::string& passName);

        /**
         * Class destructor.
         */
        virtual ~EventImpl();

    /** ********* Implementation of base class methods  ********** **/

        /**
         * Get the event header
         */
        virtual const event::EventHeader* getEventHeader() const { return eventHeader_; }

        /**
         * Adds a clones array to the event/tree.
         * @param collectionName
         * @param tca The clones array to add.
         */
        virtual void add(const std::string& collectionName, TClonesArray* tca);


       /**
         * Adds a general object to the event/tree.
         * @param name The name of the object.
         * @param obj The object to add.
         *
         * @note
         * All objects must implement/replace TObject::Clone() to
         * simply call "new" and create an empty new object and
         * implement TObject::Copy() to either copy the contents of
         * the object or swap them to the calling function, which
         * is more efficient.
         */
        virtual void add(const std::string& name, TObject* obj);

  protected:
        /**
         * Get an object from the event using a custom pass name.
         * @param collectionName The collection name.
         * @param passName The pass name.
         */
    virtual const TObject* getReal(const std::string& collectionName, const std::string& passName, bool mustExist);

  public:
    
    /** ********* Functionality for storage  ********** **/


         event::EventHeader& getEventHeaderMutable() const { return *eventHeader_; }
    
        /**
         * Set the input data tree.
         * @param tree The input data tree.
         */
        void setInputTree(TTree* tree);

        /**
         * Set the output data tree.
         * @param tree The output data tree.
         */
        void setOutputTree(TTree* tree);

        /**
         * Create the output data tree.
         * @return The output data tree.
         */
        TTree* createTree();

        /**
         * Make a branch name from a collection and pass name.
         * @param collectionName The collection name.
         * @param passName The pass name.
         */
        std::string makeBranchName(const std::string& collectionName, const std::string& passName) const {
            return collectionName + "_" + passName;
        }

        /**
         * Make a branch name from a collection and the default(current) pass name.
         * @param collectionName The collection name.
         */
        std::string makeBranchName(const std::string& collectionName) const {
            return makeBranchName(collectionName, passName_);
        }


        /*
         * These two methods only apply to the current pass of processing --
         * it is not allowed to modify any object from a previous pass.  They will
         * require null if there no object with that name in the current pass.
         */
        // TClonesArray* getMutable(const std::string& collectionName,const std::string& passName);
        // TObject* getMutable(const std::string& collectionName, const std::string& passName);


        /**
         * Go to the next event by incrementing the entry index.
         * @return Hard-coded to return true.
         */
        bool nextEvent();

        void beforeFill();
        void Clear();
    
        /**
         * Perform end of event action (clears the owned objects).
         */
        void onEndOfEvent();

        /**
         * Perform end of file action (doesn't do anything right now).
         */
        void onEndOfFile();

        /**
         * Get the current/default pass name.
         * @return The current/default pass name.
         */
        std::string getPassName() {
            return passName_;
        }

    private:

        /**
         * The event header object (as pointer)
         */
        event::EventHeader* eventHeader_{nullptr};
    
        /**
         * Number of entries in the tree.
         */
        Long64_t entries_{-1};

        /**
         * Current entry in the tree.
         */
        Long64_t ientry_ {-1};

        /**
         * The default pass name.
         */
        std::string passName_;

        /**
         * The output tree for writing a new file.
         */
        TTree* outputTree_{nullptr};

        /**
         * The input tree for reading existing data.
         */
        TTree* inputTree_{nullptr};

        /**
         * Map of names to branches.
         */
        mutable std::map<std::string, TBranch*> branches_;

        /**
         * Map of names to objects.
         */
        mutable std::map<std::string, TObject*> objects_;

        /**
         * Map of owned objects that should eventually be cleared at end of event
         * and deleted when this object is destroyed.
         */
        std::map<std::string, TObject*> objectsOwned_;

        /**
         * List of new branches added.
         */
        std::vector<TBranch*> newBranches_;


        /**
         * Names of all branches
         */
        std::vector<std::string> branchNames_;


    
        /**
         * Efficiency cache for empty pass name lookups
         */
        mutable std::map<std::string,std::string> knownLookups_;
};

}

#endif
