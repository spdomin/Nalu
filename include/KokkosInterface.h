/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef INCLUDE_KOKKOSINTERFACE_H_
#define INCLUDE_KOKKOSINTERFACE_H_

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

using HostSpace = Kokkos::HostSpace;
using DeviceSpace = Kokkos::DefaultExecutionSpace;
using DeviceShmem = DeviceSpace::scratch_memory_space;
typedef Kokkos::Schedule<Kokkos::Dynamic> DynamicScheduleType;
typedef typename Kokkos::TeamPolicy<typename Kokkos::DefaultExecutionSpace, DynamicScheduleType>::member_type TeamHandleType;
template <typename T>
using SharedMemView = Kokkos::View<T, Kokkos::LayoutRight, DeviceShmem, Kokkos::MemoryUnmanaged>;

//!!!! Important: the 'Kokkos::Serial' here makes sure we only run serial. !!!!
using DeviceTeamPolicy = Kokkos::TeamPolicy<Kokkos::Serial>;
using DeviceTeam = DeviceTeamPolicy::member_type;

inline DeviceTeamPolicy get_team_policy(const size_t sz, const size_t bytes_per_team,
    const size_t bytes_per_thread)
{
  DeviceTeamPolicy policy(sz, Kokkos::AUTO);
  return policy.set_scratch_size(0, Kokkos::PerTeam(bytes_per_team), Kokkos::PerThread(bytes_per_thread));
}

inline
SharedMemView<stk::mesh::Entity*> get_entity_shmem_view_1D(const TeamHandleType& team, size_t len)
{
  return Kokkos::subview(SharedMemView<stk::mesh::Entity**>(team.team_shmem(), team.team_size(), len), team.team_rank(), Kokkos::ALL());
}

inline
SharedMemView<double*> get_shmem_view_1D(const TeamHandleType& team, size_t len)
{
  return Kokkos::subview(SharedMemView<double**>(team.team_shmem(), team.team_size(), len), team.team_rank(), Kokkos::ALL());
}

inline
SharedMemView<double**> get_shmem_view_2D(const TeamHandleType& team, size_t len1, size_t len2)
{
  return Kokkos::subview(SharedMemView<double***>(team.team_shmem(), team.team_size(), len1, len2), team.team_rank(), Kokkos::ALL(), Kokkos::ALL());
}

inline
SharedMemView<double***> get_shmem_view_3D(const TeamHandleType& team, size_t len1, size_t len2, size_t len3)
{
  return Kokkos::subview(SharedMemView<double****>(team.team_shmem(), team.team_size(), len1, len2, len3), team.team_rank(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
}

template<typename SizeType, class Function>
void kokkos_parallel_for(const std::string& debuggingName, SizeType n, Function loop_body)
{
    Kokkos::parallel_for(debuggingName, Kokkos::RangePolicy<Kokkos::Serial>(0, n), loop_body);
}

template<typename SizeType, class Function, typename ReduceType>
void kokkos_parallel_reduce(SizeType n, Function loop_body, ReduceType& reduce, const std::string& debuggingName)
{
    Kokkos::parallel_reduce(debuggingName, Kokkos::RangePolicy<Kokkos::Serial>(0, n), loop_body, reduce);
}

}
}

#endif /* INCLUDE_KOKKOSINTERFACE_H_ */
