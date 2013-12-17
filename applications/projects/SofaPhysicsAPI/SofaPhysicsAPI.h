#ifndef SOFAPHYSICSAPI_H
#define SOFAPHYSICSAPI_H

class SofaPhysicsSimulation;
class SofaPhysicsOutputMesh;
class SofaPhysicsOutputMeshTetrahedron;
class SofaPhysicsDataMonitor;
class SofaPhysicsDataController;

typedef unsigned int Index; ///< Type used for topology indices
typedef float Real;         ///< Type used for coordinates
typedef void* ID;           ///< Type used for IDs

/// Main class used to control a Sofa Simulation
class SofaPhysicsSimulation
{
public:
    SofaPhysicsSimulation(bool useGUI = false, int GUIFramerate = 0);
    ~SofaPhysicsSimulation();

    /// Load an XML file containing the main scene description
    bool load(const char* filename);

    /// Start the simulation
    /// Currently this simply sets the animated flag to true, but this might
    /// start a separate computation thread in a future version
    void start();

    /// Stop/pause the simulation
    void stop();

    /// Compute one simulation time-step
    void step();

    /// Reset the simulation to its initial state
    void reset();

    /// Send an event to the simulation for custom controls
    /// (such as switching active instrument)
    void sendValue(const char* name, double value);

    /// Reset the camera to its default position
    void resetView();

    /// Render the scene using OpenGL
    void drawGL();

    /// Return the number of currently active output meshes
    unsigned int            getNbOutputMeshes();

    /// Return the number of currently active output Tetrahedron meshes
    unsigned int getNbOutputMeshTetrahedrons();

    /// Return an array of pointers to active output meshes
    SofaPhysicsOutputMesh** getOutputMeshes();

    SofaPhysicsOutputMeshTetrahedron** getOutputMeshTetrahedrons();

    /// Return true if the simulation is running
    /// Note that currently you must call the step() method
    /// periodically to actually animate the scene
    bool isAnimated() const;

    /// Set the animated state to a given value (requires a
    /// simulation to be loaded)
    void setAnimated(bool val);

    /// Return the main simulation file name (from the last
    /// call to load())
    const char* getSceneFileName() const;

    /// Return the current time-step (or 0 if no simulation
    /// is loaded)
    double getTimeStep() const;
    /// Control the timestep of the simulation (requires a
    /// simulation to be loaded)
    void   setTimeStep(double dt);

    /// Return the current simulated time
    double getTime() const;

    /// Return the current computation speed (averaged over
    /// the last 100 steps)
    double getCurrentFPS() const;

    /// Return the number of currently active data monitors
    unsigned int getNbDataMonitors();

    /// Return an array of pointers to active data monitors
    SofaPhysicsDataMonitor** getDataMonitors();

    /// Return the number of currently active data controllers
    unsigned int getNbDataControllers();

    /// Return an array of pointers to active data controllers
    SofaPhysicsDataController** getDataControllers();

    /// Internal implementation sub-class
    class Impl;
    /// Internal implementation sub-class
    Impl* impl;
};

/// Class describing one output mesh (i.e. visual model) in the simulation
class SofaPhysicsOutputMesh
{
public:

    SofaPhysicsOutputMesh();
    ~SofaPhysicsOutputMesh();

    const char* getName(); ///< (non-unique) name of this object
    ID          getID();   ///< unique ID of this object

    unsigned int getNbVertices(); ///< number of vertices
    const Real* getVPositions();  ///< vertices positions (Vec3)
    const Real* getVNormals();    ///< vertices normals   (Vec3)
    const Real* getVTexCoords();  ///< vertices UVs       (Vec2)
    int getTexCoordRevision();    ///< changes each time texture coord data are updated
    int getVerticesRevision();    ///< changes each time vertices data are updated

    unsigned int getNbVAttributes();                    ///< number of vertices attributes
    unsigned int getNbAttributes(int index);            ///< number of the attributes in specified vertex attribute 
    const char*  getVAttributeName(int index);          ///< vertices attribute name
    int          getVAttributeSizePerVertex(int index); ///< vertices attribute #
    const Real*  getVAttributeValue(int index);         ///< vertices attribute (Vec#)
    int          getVAttributeRevision(int index);      ///< changes each time vertices attribute is updated

    unsigned int getNbLines(); ///< number of lines
    const Index* getLines();   ///< lines topology (2 indices / line)
    int getLinesRevision();    ///< changes each time lines data is updated

    unsigned int getNbTriangles(); ///< number of triangles
    const Index* getTriangles();   ///< triangles topology (3 indices / triangle)
    int getTrianglesRevision();    ///< changes each time triangles data is updated

    unsigned int getNbQuads(); ///< number of quads
    const Index* getQuads();   ///< quads topology (4 indices / quad)
    int getQuadsRevision();    ///< changes each time quads data is updated

    /// Internal implementation sub-class
    class Impl;
    /// Internal implementation sub-class
    Impl* impl;
};

class SofaPhysicsOutputMeshTetrahedron 
{
public:

  SofaPhysicsOutputMeshTetrahedron();
  ~SofaPhysicsOutputMeshTetrahedron();

  const char* getName(); ///< (non-unique) name of this object
  ID          getID();   ///< unique ID of this object

  unsigned int getNbVertices(); ///< number of vertices
  const Real* getVPositions();  ///< vertices positions (Vec3)
  const Real* getVNormals();    ///< vertices normals   (Vec3)
  const Real* getVTexCoords();  ///< vertices UVs       (Vec2)
  int getTexCoordRevision();    ///< changes each time texture coord data are updated
  int getVerticesRevision();    ///< changes each time vertices data are updated

  unsigned int getNbVAttributes();                    ///< number of vertices attributes
  unsigned int getNbAttributes(int index);            ///< number of the attributes in specified vertex attribute 
  const char*  getVAttributeName(int index);          ///< vertices attribute name
  int          getVAttributeSizePerVertex(int index); ///< vertices attribute #
  const Real*  getVAttributeValue(int index);         ///< vertices attribute (Vec#)
  int          getVAttributeRevision(int index);      ///< changes each time vertices attribute is updated

  unsigned int getNbLines(); ///< number of lines
  const Index* getLines();   ///< lines topology (2 indices / line)
  int getLinesRevision();    ///< changes each time lines data is updated

  unsigned int getNbTriangles(); ///< number of triangles
  const Index* getTriangles();   ///< triangles topology (3 indices / triangle)
  int getTrianglesRevision();    ///< changes each time triangles data is updated

  unsigned int getNbQuads(); ///< number of quads
  const Index* getQuads();   ///< quads topology (4 indices / quad)
  int getQuadsRevision();    ///< changes each time quads data is updated

  unsigned int getNbTetrahedrons(); ///< number of Tetrahedrons
  const Index* getTetrahedrons();   ///< Tetrahedrons topology (4 indices / Tetrahedron)
  int getTetrahedronsRevision();    ///< changes each time Tetrahedrons data is updated

  /// Internal implementation sub-class
  class Impl;
  /// Internal implementation sub-class
  Impl* impl;
};

/// Class for data monitoring
class SofaPhysicsDataMonitor
{
public:

    SofaPhysicsDataMonitor();
    ~SofaPhysicsDataMonitor();

    const char* getName(); ///< (non-unique) name of this object
    ID          getID();   ///< unique ID of this object

    /// Get the value of the associated variable
    const char* getValue();

    /// Internal implementation sub-class
    class Impl;
    /// Internal implementation sub-class
    Impl* impl;
};

/// Class for data control
class SofaPhysicsDataController
{
public:

    SofaPhysicsDataController();
    ~SofaPhysicsDataController();

    const char* getName(); ///< (non-unique) name of this object
    ID          getID();   ///< unique ID of this object

    /// Set the value of the associated variable
    void setValue(const char* v);

    /// Internal implementation sub-class
    class Impl;
    /// Internal implementation sub-class
    Impl* impl;
};

#endif // SOFAPHYSICSAPI_H
