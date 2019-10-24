//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_INTERNALS_CARTESIAN_RAY_H_INCLUDED )
#define  KRATOS_INTERNALS_CARTESIAN_RAY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{
namespace Internals
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// A cartesian ray to be used in ray casting operations
/** This class represents a cartesian ray in 3D space defined by a direction and Point1 and Point2
*/
template<typename TGeometryType> 
class CartesianRay
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CartesianRay
    KRATOS_CLASS_POINTER_DEFINITION(CartesianRay);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CartesianRay(): mDirection(0), mPoint1(), mPoint2() {}

    /// Constructor with all needed parameters
    CartesianRay(int Direction, Point const& Point1, Point const& Point2): mDirection(Direction), mPoint1(Point1), mPoint2(Point2) {}

    // Copy constructor
    CartesianRay(CartesianRay const& Other): mDirection(Other.mDirection), mPoint1(Other.mPoint1), mPoint2(Other.mPoint2), mIntersections(Other.mIntersections) {}

    /// Destructor.
    virtual ~CartesianRay(){}

    ///@}
    ///@name Operators
    ///@{

    CartesianRay& operator=(CartesianRay const& Other){
        mDirection=Other.mDirection;
        mPoint1 = Other.mPoint1;
        mPoint2 = Other.mPoint2;
        mIntersections = Other.mIntersections;

        return *this;
    }


    ///@}
    ///@name Operations
    ///@{

    void AddIntersection(TGeometryType const& rGeometry, double Tolerance){

        array_1d<double,3> intersection_point = ZeroVector(3);
        std::vector<double> geometry_point_1{rGeometry[0].X(), rGeometry[0].Y(), rGeometry[0].Z()};
        std::vector<double> geometry_point_2{rGeometry[1].X(), rGeometry[1].Y(), rGeometry[1].Z()};
        std::vector<double> geometry_point_3{rGeometry[2].X(), rGeometry[2].Y(), rGeometry[2].Z()};
        
        std::vector<double> point_1{mPoint1.X(), mPoint1.Y(), mPoint1.Z()};
        std::vector<double> point_2{mPoint2.X(), mPoint2.Y(), mPoint2.Z()};

        Point min_point;
        Point max_point;
        max_point = *(rGeometry.begin());
        min_point = *(rGeometry.begin());
        for(auto const& point : rGeometry){
            for(std::size_t i = 0; i<3; i++)
            {
                min_point[i] =  (min_point[i] >  point[i] ) ?  point[i] : min_point[i];
                max_point[i] =  (max_point[i] <  point[i] ) ?  point[i] : max_point[i];
            }
        }

		const double relative_tolerance = 1.0e-12*std::sqrt(rGeometry.Length());
        const int is_intersected = IntersectionUtilities::ComputeTriangleLineIntersection(
          rGeometry,
          mPoint1,
          mPoint2,
          intersection_point,
          relative_tolerance);

        if(is_intersected == 1){ // There is an intersection but not coplanar
            mIntersections.push_back(std::make_pair(intersection_point[mDirection], &rGeometry));
        }
    }

    void CollapseIntersectionPoints(double Tolerance){

        if (mIntersections.empty()) {
            return;
        }

        std::size_t new_size = 0;
        // Sort
        std::sort(mIntersections.begin(), mIntersections.end());
        // Unique
        auto i_begin = mIntersections.begin();
        auto i_intersection = mIntersections.begin();
        while (++i_begin != mIntersections.end()) {
            // considering the very near points as the same points
            if (std::abs(i_begin->first - i_intersection->first) > Tolerance) { // if the hit points are far enough they are not the same 
                if(new_size == 0)
                    new_size++; // we should consider the first intersection then.               
                *(++i_intersection) = *i_begin;
                new_size++;
            }
            else{ // Now there are near hits, so we check if it is really pass through a duplicated surface or just passing tangent to model
                // Getting the patch of geometries with near hit
                auto i_patch_begin = i_intersection;
                auto i_patch_end = i_begin;
                while(++i_begin != mIntersections.end()){
                    if (std::abs(i_begin->first - i_intersection->first) > Tolerance) { // This is the end of the patch.               
                        break;
                    }
                }
                
                i_patch_end = i_begin;

                if(CheckPassingThroughByExtraRays(i_patch_begin, i_patch_end, Tolerance,  2.00*Tolerance)) {
                    *(++i_intersection) = *i_begin;
                    new_size++;
                }

                if(i_begin == mIntersections.end())
                    break;
            }
        }
        mIntersections.resize(new_size);
    }

    void CalculateColor(std::vector<double> const& Coordinates, int InsideColor, int OutsideColor, std::vector<double>& ResultingColors, double NearEnough){

        bool is_inside=false;

        if(ResultingColors.size() != Coordinates.size())
        ResultingColors.resize(Coordinates.size());
    
        std::size_t current_index=0;
        const std::size_t size = Coordinates.size();

        for(auto& i_intersection : mIntersections){
            while(current_index < size){
                if((i_intersection.first - Coordinates[current_index]) > NearEnough){
                    ResultingColors[current_index++] = (is_inside) ? InsideColor : OutsideColor;
                } else if((i_intersection.first - Coordinates[current_index]) > -NearEnough){ // Considering very near to wall as inside
                    ResultingColors[current_index++] =  InsideColor;
                }
                else{
                    break;
                }
            }
            is_inside = (is_inside) ? false:true;
        }
        // now coloring the points after the last intersection as outside
        while(current_index < size){
            ResultingColors[current_index++] = OutsideColor;
        }
    }


    ///@}
    ///@name Access
    ///@{

    std::vector<std::pair<double, const TGeometryType*>> const& GetIntersections() const {return mIntersections;}

    Point& GetPoint1(){return mPoint1;}

    Point& GetPoint2(){return mPoint2;}


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "CartesianRay" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "CartesianRay";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    int mDirection;
    Point mPoint1;
    Point mPoint2;
    std::vector<std::pair<double, const TGeometryType*>> mIntersections;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    bool CheckPassingThroughByExtraRays(typename std::vector<std::pair<double, const TGeometryType*>>::iterator Begin, typename std::vector<std::pair<double, const TGeometryType*>>::iterator End, double Tolerance, double Delta){
        std::array<double, 8> delta_u{Delta, Delta, 0.00, -Delta, -Delta, -Delta, 0.00, Delta};
        std::array<double, 8> delta_v{0.00, Delta, Delta, Delta, 0.00, -Delta, -Delta, -Delta};

        std::size_t no_hit_cases = 0;

        for(int i_ray = 0 ; i_ray < 8 ; i_ray++){
            CartesianRay extra_ray(mDirection, mPoint1, mPoint2);
            if(mDirection == 0) {
                extra_ray.mPoint1[1] += delta_u[i_ray];
                extra_ray.mPoint1[2] += delta_v[i_ray];
                extra_ray.mPoint2[1] += delta_u[i_ray];
                extra_ray.mPoint2[2] += delta_v[i_ray];
            }
            else if(mDirection == 1) {
                extra_ray.mPoint1[0] += delta_u[i_ray];
                extra_ray.mPoint1[2] += delta_v[i_ray];
                extra_ray.mPoint2[0] += delta_u[i_ray];
                extra_ray.mPoint2[2] += delta_v[i_ray];
            }
            else if(mDirection == 2) {
                extra_ray.mPoint1[0] += delta_u[i_ray];
                extra_ray.mPoint1[1] += delta_v[i_ray];
                extra_ray.mPoint2[0] += delta_u[i_ray];
                extra_ray.mPoint2[1] += delta_v[i_ray];
            }
            for(auto i_intersection = Begin ; i_intersection != End; i_intersection++){
                extra_ray.AddIntersection(*(i_intersection->second), Tolerance);
            }
            if(extra_ray.mIntersections.size() == 0){
                no_hit_cases++;
            }
            if(no_hit_cases > 4) // more than half
                return false;
        }
        return true;
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class CartesianRay

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<typename TGeometryType> 
inline std::istream& operator >> (std::istream& rIStream,
                CartesianRay<TGeometryType>& rThis){
                    return rIStream;
                }

/// output stream function
template<typename TGeometryType> 
inline std::ostream& operator << (std::ostream& rOStream,
                const CartesianRay<TGeometryType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Internals

}  // namespace Kratos.

#endif // KRATOS_INTERNALS_CARTESIAN_RAY_H_INCLUDED  defined

