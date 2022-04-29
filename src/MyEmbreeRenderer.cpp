
// MODIFIED VERSION OF LIBIGL'S EmbreeRenderer.cpp

#include "MyEmbreeRenderer.h"
// Implementation

//IGL viewing parts
#include <igl/unproject.h>
#include <igl/look_at.h>
#include <igl/frustum.h>
#include <igl/ortho.h>

// color map
#include <igl/jet.h>

#include <igl/PI.h>


IGL_INLINE void igl::embree::EmbreeRenderer::init_view()
{
  camera_base_zoom = 1.0f;
  camera_zoom = 1.0f;

  camera_view_angle = 45.0;
  camera_dnear = 1.0;
  camera_dfar = 100.0;
  camera_base_translation << 0, 0, 0;
  camera_translation << 0, 0, 0;
  camera_eye << 0, 0, 5;
  camera_center << 0, 0, 0;
  camera_up << 0, 1, 0;
  
  rot_matrix = Eigen::Matrix3f::Identity();

  view = Eigen::Matrix4f::Identity();
  proj = Eigen::Matrix4f::Identity();
  norm = Eigen::Matrix4f::Identity();

  orthographic = false;
  

  uC << 1,0,0;
}

IGL_INLINE igl::embree::EmbreeRenderer::EmbreeRenderer()
  :
  scene(NULL),
  geomID(0),
  initialized(false),
  device(igl::embree::EmbreeDevice::get_device())
{
  init_view();
}

IGL_INLINE igl::embree::EmbreeRenderer::EmbreeRenderer(
  const EmbreeRenderer &)
  :// To make -Weffc++ happy
  scene(NULL),
  geomID(0),
  initialized(false)
{
  assert(false && "Embree: Copying EmbreeRenderer is not allowed");
}

IGL_INLINE igl::embree::EmbreeRenderer & igl::embree::EmbreeRenderer::operator=(
  const EmbreeRenderer &)
{
  assert(false && "Embree: Assigning an EmbreeRenderer is not allowed");
  return *this;
}


IGL_INLINE void igl::embree::EmbreeRenderer::init(
  const PointMatrixType& V,
  const FaceMatrixType& F,
  bool isStatic)
{
  std::vector<const PointMatrixType*> Vtemp;
  std::vector<const FaceMatrixType*> Ftemp;
  std::vector<int> masks;
  Vtemp.push_back(&V);
  Ftemp.push_back(&F);
  masks.push_back(0xFFFFFFFF);
  init(Vtemp,Ftemp,masks,isStatic);
}

IGL_INLINE void igl::embree::EmbreeRenderer::init(
  const std::vector<const PointMatrixType*>& V,
  const std::vector<const FaceMatrixType*>&  F,
  const std::vector<int>& masks,
  bool isStatic)
{
  if(initialized)
    deinit();

  using namespace std;

  if(V.size() == 0 || F.size() == 0)
  {
    std::cerr << "Embree: No geometry specified!";
    return;
  }
  RTCBuildQuality buildQuality = isStatic ? RTC_BUILD_QUALITY_HIGH : RTC_BUILD_QUALITY_MEDIUM;

  // create a scene
  scene = rtcNewScene(device);
  rtcSetSceneFlags(scene, RTC_SCENE_FLAG_ROBUST);
  rtcSetSceneBuildQuality(scene, buildQuality);

  for(int g=0;g<(int)V.size();g++)
  {
    // create triangle mesh geometry in that scene
    RTCGeometry geom_0 = rtcNewGeometry (device, RTC_GEOMETRY_TYPE_TRIANGLE);
    rtcSetGeometryBuildQuality(geom_0, buildQuality);
    rtcSetGeometryTimeStepCount(geom_0,1);
    geomID = rtcAttachGeometry(scene,geom_0);
    rtcReleaseGeometry(geom_0);

    // fill vertex buffer, have to be 16 byte wide( sizeof(float)*4 )
    Eigen::Map<Eigen::Matrix<float,-1,4,Eigen::RowMajor>> vertices(
        (float*)rtcSetNewGeometryBuffer(geom_0,RTC_BUFFER_TYPE_VERTEX,0,RTC_FORMAT_FLOAT3,4*sizeof(float),V[g]->rows()),
        V[g]->rows(),4
    );
    vertices.block(0,0,V[g]->rows(),3) = V[g]->cast<float>();

    // fill triangle buffer
    Eigen::Map<Eigen::Matrix<unsigned int,-1,3,Eigen::RowMajor>> triangles(
      (unsigned int*) rtcSetNewGeometryBuffer(geom_0,RTC_BUFFER_TYPE_INDEX,0,RTC_FORMAT_UINT3,3*sizeof(unsigned int), F[g]->rows()),
      F[g]->rows(),3
    );
    triangles = F[g]->cast<unsigned int>();
    //TODO: store vertices and triangles in array for whatever reason?

    rtcSetGeometryMask(geom_0, masks[g]);
    rtcCommitGeometry(geom_0);
  }

  rtcCommitScene(scene);

  if(rtcGetDeviceError (device) != RTC_ERROR_NONE)
      std::cerr << "Embree: An error occurred while initializing the provided geometry!" << endl;
#ifdef IGL_VERBOSE
  else
    std::cerr << "Embree: geometry added." << endl;
#endif

  initialized = true;
}

IGL_INLINE igl::embree::EmbreeRenderer::~EmbreeRenderer()
{
  if(initialized)
    deinit();
  
  igl::embree::EmbreeDevice::release_device();
}

IGL_INLINE void igl::embree::EmbreeRenderer::deinit()
{
  if(scene)
  {
    rtcReleaseScene(scene);

    if(rtcGetDeviceError (device) != RTC_ERROR_NONE)
    {
        std::cerr << "Embree: An error occurred while resetting!" << std::endl;
    }
#ifdef IGL_VERBOSE
    else
    {
      std::cerr << "Embree: geometry removed." << std::endl;
    }
#endif
  }
}

IGL_INLINE bool igl::embree::EmbreeRenderer::intersect_ray(
  const Eigen::RowVector3f& origin,
  const Eigen::RowVector3f& direction,
  Hit&  hit,
  float tnear,
  float tfar,
  int mask) const
{
  RTCRayHit ray;
  ray.ray.flags = 0;
  create_ray(ray, origin,direction,tnear,tfar,mask);

  // shot ray
  {
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);
    rtcIntersect1(scene, &context, &ray);
    ray.hit.Ng_x = -ray.hit.Ng_x; // EMBREE_FIXME: only correct for triangles,quads, and subdivision surfaces
    ray.hit.Ng_y = -ray.hit.Ng_y;
    ray.hit.Ng_z = -ray.hit.Ng_z;
  }
#ifdef IGL_VERBOSE
  if(rtcGetDeviceError (device) != RTC_ERROR_NONE)
      std::cerr << "Embree: An error occurred while resetting!" << std::endl;
#endif

  if((unsigned)ray.hit.geomID != RTC_INVALID_GEOMETRY_ID)
  {
    hit.id = ray.hit.primID;
    hit.gid = ray.hit.geomID;
    hit.u = ray.hit.u;
    hit.v = ray.hit.v;
    hit.t = ray.ray.tfar;
    hit.N = Vec3f(ray.hit.Ng_x, ray.hit.Ng_y, ray.hit.Ng_z);
    return true;
  }

  return false;
}

IGL_INLINE void 
igl::embree::EmbreeRenderer
::create_ray(RTCRayHit& ray, const Eigen::RowVector3f& origin, 
 const Eigen::RowVector3f& direction, float tnear, float tfar, int mask) const
{
  ray.ray.org_x = origin[0];
  ray.ray.org_y = origin[1];
  ray.ray.org_z = origin[2];
  ray.ray.dir_x = direction[0];
  ray.ray.dir_y = direction[1];
  ray.ray.dir_z = direction[2];
  ray.ray.tnear = tnear;
  ray.ray.tfar = tfar;
  ray.ray.id = RTC_INVALID_GEOMETRY_ID;
  ray.ray.mask = mask;
  ray.ray.time = 0.0f;

  ray.hit.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
  ray.hit.primID = RTC_INVALID_GEOMETRY_ID;
}


IGL_INLINE void
igl::embree::EmbreeRenderer
::set_mesh(const Eigen::Matrix<double,Eigen::Dynamic,3> & MV,
           const Eigen::Matrix<int,  Eigen::Dynamic,3>  & MF,
           bool is_static)
{
  V = MV.cast<float>();
  F = MF;
  this->init(V,F,is_static);

  auto min_point = V.colwise().minCoeff();
  auto max_point = V.colwise().maxCoeff();
  auto centroid  = (0.5*(min_point + max_point)).eval();

  camera_base_translation.setConstant(0);
  camera_base_translation.head(centroid.size()) = -centroid.cast<float>();
  camera_base_zoom = 2.5 / (max_point-min_point).array().abs().maxCoeff();
}

IGL_INLINE void
igl::embree::EmbreeRenderer
::render_buffer(PixelMatrixType& R, PixelMatrixType&G, PixelMatrixType &B,PixelMatrixType &A)
{
  assert(R.rows()==G.rows());assert(R.rows()==B.rows());assert(R.rows()==A.rows());
  assert(R.cols()==G.cols());assert(R.cols()==B.cols());assert(R.cols()==A.cols());

  Eigen::Vector4f viewport(0,0,R.rows(),R.cols());

  float width  = R.rows();
  float height = R.cols();

  // update view matrix
  igl::look_at( camera_eye, camera_center, camera_up, view);

  view = view
    * (rot_matrix * Eigen::Scaling(camera_zoom * camera_base_zoom)
    * Eigen::Translation3f(camera_translation + camera_base_translation)).matrix();
  
  if (orthographic)
  {
    float length = (camera_eye - camera_center).norm();
    float h = tan(camera_view_angle/360.0 * igl::PI) * (length);
    igl::ortho(-h*width/height, h*width/height, -h, h, camera_dnear, camera_dfar, proj);
  } else {
    float fH = tan(camera_view_angle / 360.0 * igl::PI) * camera_dnear;
    float fW = fH * (double)width/(double)height;
    igl::frustum(-fW, fW, -fH, fH, camera_dnear, camera_dfar, proj);
  }

  // go over all pixels in the "view"
  for(int x=0;x<(int)width;++x) 
  {
    for(int y=0;y<(int)height;++y) 
    {
      Vec3f s,d,dir;
      igl::embree::EmbreeRenderer::Hit hit;

      // look into the screen
      Vec3f win_s(x,y,0);
      Vec3f win_d(x,y,1);
      // Source, destination and direction in world
      igl::unproject(win_s,this->view,this->proj,viewport,s);
      igl::unproject(win_d,this->view,this->proj,viewport,d);
      dir = d-s;
      dir.normalize();

      auto clamp=[](float x)->unsigned char {return (unsigned char)(x<0?0:x>1.0?255:x*255);};

      if(this->intersect_ray(s,dir,hit))
      {
        if (true || dir.dot(hit.N) > 0.0f) // "true": no front face culling
        {
          if (dir.dot(hit.N) < 0.0f) {
            hit.N = - hit.N;
          }
          // TODO: interpolate normals ?
          hit.N.normalize();
          // cos between ray and face normal
          float face_proj=dir.dot(hit.N);
          face_proj = 0.4 + 1.0*face_proj/2.0;
          face_proj *= 1.5; // TODO Parameter with this
          //face_proj = std::sqrt(face_proj);

          face_proj -= hit.t/1000.0;
          //face_proj = std::sqrt(face_proj);

          Eigen::RowVector3f c;

          if(this->uniform_color)
          {
            // same color for the whole mesh
            c=uC;
          } else if(this->face_based) {
            // flat color per face
            c=this->C.row(hit.id);
          } else { //use barycentric coordinates to interpolate colour
              c=this->C.row(F(hit.id,1))*hit.u+
                this->C.row(F(hit.id,2))*hit.v+
                this->C.row(F(hit.id,0))*(1.0-hit.u-hit.v);
          }

          R(x,y) = clamp(face_proj*c(0));
          G(x,y) = clamp(face_proj*c(1));
          B(x,y) = clamp(face_proj*c(2));

          double threshold = 0.05;
          Eigen::RowVector3f p = V.row(F(hit.id, 1)) * hit.u 
                               + V.row(F(hit.id, 2)) * hit.v 
                               + V.row(F(hit.id, 0)) * (1 - hit.u - hit.v); // WHY THE F ARE BARYCENTRIC DONE SUCH THAT u CORRESPONDS TO 1 AND NOT 0 ?????!!!!

          double best_dist = -1;

          if (hit.id < edges_per_tri_.size()) {
            for (int edge_id: edges_per_tri_[hit.id]){
              
              // p point
              // p1 p2 edge of interest
              // pp projection of p on p1p2
              
              Eigen::RowVector3f p1 = begs_.row(edge_id);
              Eigen::RowVector3f p2 = ends_.row(edge_id);
              
              float dotted = (p - p1).dot( (p2 - p1).normalized() );
              Eigen::RowVector3f ppp1 = dotted * (p2 - p1) / (p2 - p1).norm();
              Eigen::RowVector3f pp = p1 + ppp1;

              //if (ppp1.norm() < 0 || ppp1.norm() > (p2 - p1).norm()) continue;

              double dist = (p - pp).norm();

              if (dotted < 0){
                dist = (p - p1).norm();
              }
              else if (dotted > (p2 - p1).norm()){
                dist = (p - p2).norm();
              }

              if (dist < threshold && (best_dist == -1 || dist < best_dist)){
                R(x,y) = 255 * edge_colors_(edge_id, 0);
                G(x,y) = 255 * edge_colors_(edge_id, 1);
                B(x,y) = 255 * edge_colors_(edge_id, 2);
                A(x,y) = 255;
                best_dist = dist;
              }
            }
          }


          /*
          for (int e=0; e<3; e++){
            int v1 = F(hit.id, e);
            if (edge_colors_.row(v1)(0) == -1) continue;
            for (int e2=e+1; e2<3; e2++){
              int v2 = F(hit.id, e2);
              if (edge_colors_.row(v2)(0) == -1) continue;

              Eigen::RowVector3f A = V.row(v1);
              Eigen::RowVector3f B = V.row(F(hit.id, (e2+1)%3));
              // from p to edge
              double dist = ((A - p) + (B - A) * (p - A).dot(B - A) / ((B - A).norm() * (B - A).norm())).norm(); 

              if (dist < threshold){
                R(x,y) = 0;
                G(x,y) = 0;
                B(x,y) = 0;
                A(x,y) = 0;
              }
            }
            
          }*/
        }
        // give the same alpha to all points with something behind
        A(x,y)=255; 
      } else {
        R(x,y)=255;
        G(x,y)=255;
        B(x,y)=255;
        A(x,y)=0; // transparent background here
      }

    }
  }
}


IGL_INLINE void
igl::embree::EmbreeRenderer
::set_colors(const Eigen::MatrixXd & C)
{
  if(C.rows()==V.rows()) // per vertex color
  {
    face_based = false;
    this->C = C.cast<float>();
    this->uniform_color=false;
  } else if (C.rows()==F.rows()) {
    face_based = true;
    this->C = C.cast<float>();
    this->uniform_color=false;
  } else if (C.rows()==1) {
    face_based = true;
    this->uC = C.cast<float>();
    this->uniform_color=true;
  }else {
    // don't know what to do
    this->uniform_color=true;
    assert(false); //?
  }
}

IGL_INLINE void
igl::embree::EmbreeRenderer
::set_data(const Eigen::VectorXd & D, igl::ColorMapType cmap)
{
  const double caxis_min = D.minCoeff();
  const double caxis_max = D.maxCoeff();
  return set_data(D,caxis_min,caxis_max,cmap);
}


IGL_INLINE void igl::embree::EmbreeRenderer::set_data(
  const Eigen::VectorXd & D,
  double caxis_min,
  double caxis_max,
  igl::ColorMapType cmap)
{
    Eigen::MatrixXd C;
    igl::colormap(cmap,D,caxis_min,caxis_max,C);
    set_colors(C);
}

IGL_INLINE void 
igl::embree::EmbreeRenderer::set_rot(const Eigen::Matrix3d &r)
{
  this->rot_matrix = r.cast<float>();
}

IGL_INLINE void 
igl::embree::EmbreeRenderer::set_zoom(double zoom)
{
  this->camera_zoom=zoom;
}

IGL_INLINE void 
igl::embree::EmbreeRenderer::set_translation(const Eigen::Vector3d &tr)
{
  this->camera_translation=tr.cast<float>();
}

IGL_INLINE void 
igl::embree::EmbreeRenderer::set_face_based(bool _f)
{
  this->face_based=_f;
}

IGL_INLINE void 
igl::embree::EmbreeRenderer::set_orthographic(bool o)
{
  this->orthographic=o;
}

#ifdef IGL_STATIC_LIBRARY
#endif //IGL_STATIC_LIBRARY
