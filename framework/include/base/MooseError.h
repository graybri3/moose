/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef MOOSEERROR_H
#define MOOSEERROR_H

#include "Moose.h"
#include "MooseException.h"
#include "libmesh/threads.h"

// libMesh includes
#include "libmesh/print_trace.h"
#include "libmesh/libmesh_common.h"

// C++ includes
#include <cstdlib>

/// Application abort macro. Uses MPI_Abort if available, std::abort otherwise
#if defined(LIBMESH_HAVE_MPI)
#define MOOSE_ABORT do { MPI_Abort(libMesh::GLOBAL_COMM_WORLD, 1); std::abort(); } while (0)
#else
#define MOOSE_ABORT do { std::abort(); } while (0)
#endif

#define mooseDoOnce(do_this) do { static bool did_this_already = false; if (!did_this_already) { did_this_already = true; do_this; } } while (0)

#define mooseCheckMPIErr(err) do { if (err != MPI_SUCCESS) { if (libMesh::global_n_processors() == 1) print_trace(); libmesh_here(); MOOSE_ABORT; } } while (0)

#define mooseException(msg)                                                         \
  do                                                                                \
  {                                                                                 \
    std::ostringstream _exception_oss_;                                             \
    _exception_oss_ << msg;                                                         \
                                                                                    \
    throw MooseException(_exception_oss_.str());                                    \
  } while (0)


#ifdef NDEBUG
#define mooseAssert(asserted, msg)
#else
#define mooseAssert(asserted, msg)                                                  \
  do                                                                                \
  {                                                                                 \
    if (!(asserted))                                                                \
    {                                                                               \
      std::ostringstream _assert_oss_;                                              \
      _assert_oss_                                                                  \
        << COLOR_RED                                                                \
        << "\n\nAssertion `" #asserted "' failed\n"                                 \
        << msg                                                                      \
        << "\nat "                                                                  \
        << __FILE__ << ", line " << __LINE__                                        \
        << COLOR_DEFAULT                                                            \
        << std::endl;                                                               \
      if (Moose::_throw_on_error)                                                   \
        throw std::runtime_error(_assert_oss_.str());                               \
      else                                                                          \
      {                                                                             \
        Moose::err << _assert_oss_.str() << std::flush;                             \
        if (libMesh::global_n_processors() == 1)                                    \
          print_trace();                                                            \
        else                                                                        \
          libMesh::write_traceout();                                                \
        libmesh_here();                                                             \
        MOOSE_ABORT;                                                                \
      }                                                                             \
    }                                                                               \
  } while (0)
#endif

template <typename... Args>
[[ noreturn ]] void mooseError2(Args && ... args);

namespace moose
{
namespace internal
{

std::string mooseMsgFmt(const std::string & msg, const std::string & title, const std::string & color);

[[ noreturn ]] void mooseErrorRaw(std::string msg, const std::string prefix = "");

void mooseStreamAll(std::ostringstream & ss);

template <typename T, typename... Args>
void mooseStreamAll(std::ostringstream & ss, T && val, Args && ... args)
{
  ss << val;
  mooseStreamAll(ss, std::forward<Args>(args)...);
}

template <typename S, typename... Args>
void
mooseWarningStream(S & oss, Args && ... args)
{
  if (Moose::_warnings_are_errors)
    mooseError2(std::forward<Args>(args)...);

  std::ostringstream ss;
  mooseStreamAll(ss, args...);
  std::string msg = mooseMsgFmt(ss.str(), "*** Warning ***", COLOR_YELLOW);
  if (Moose::_throw_on_error)
    throw std::runtime_error(msg);

  oss << msg << std::flush;
}

template <typename S, typename... Args>
void
mooseInfoStream(S & oss, Args && ... args)
{
  mooseDoOnce(
    {
      std::ostringstream ss;
      mooseStreamAll(ss, args...);
      std::string msg = mooseMsgFmt(ss.str(), "*** Info ***", COLOR_CYAN);
      oss << msg << std::flush;
    }
  );
}

template <typename S, typename... Args>
void
mooseDeprecatedStream(S & oss, Args && ... args)
{
  if (Moose::_deprecated_is_error)
    mooseError2("\n\nDeprecated code:\n", std::forward<Args>(args)...);

  mooseDoOnce(
    std::ostringstream ss;
    mooseStreamAll(ss, args...);
    std::string msg = mooseMsgFmt(ss.str(), "*** Warning, This code is deprecated and will be removed in future versions!", COLOR_YELLOW);
    oss << msg;
    ss.str("");
    if (libMesh::global_n_processors() == 1)
      print_trace(ss);
    else
      libMesh::write_traceout();
    oss << ss.str() << std::endl;
  );
}

} // namespace internal
} // namespace moose

template <typename... Args>
[[ noreturn ]] void
mooseError2(Args && ... args)
{
  std::ostringstream oss;
  moose::internal::mooseStreamAll(oss, std::forward<Args>(args)...);
  moose::internal::mooseErrorRaw(oss.str());
}

template <typename... Args>
void
mooseWarning2(Args && ... args)
{
  moose::internal::mooseWarningStream(Moose::out, std::forward<Args>(args)...);
}

template <typename... Args>
void
mooseInfo2(Args && ... args)
{
  moose::internal::mooseInfoStream(Moose::out, std::forward<Args>(args)...);
}

template <typename... Args>
void
mooseDeprecated2(Args && ... args)
{
  moose::internal::mooseDeprecatedStream(Moose::out, std::forward<Args>(args)...);
}

/**
 * Global alias to Moose::out. This is here so that _console can be used anywhere
 * throughout the framework (namely by the mooseWarning() macro). If there's a local
 * object _console, that will be used, if not, we'll fall back to just using Moose::out.
 */
extern OStreamProxy & _console;

/**
 * MOOSE wrapped versions of useful libMesh macros (see libmesh_common.h)
 */
#define mooseError(msg)                                                             \
  do                                                                                \
  {                                                                                 \
    std::ostringstream _error_oss_;                                                 \
    _error_oss_ << "\n\n"                                                           \
                << COLOR_RED                                                        \
                << "\n\n*** ERROR ***\n"                                            \
                << msg                                                              \
                << COLOR_DEFAULT                                                    \
                << "\n\n";                                                          \
    if (Moose::_throw_on_error)                                                     \
      throw std::runtime_error(_error_oss_.str());                                  \
    else                                                                            \
    {                                                                               \
      Moose::err << _error_oss_.str() << std::flush;                                \
      if (libMesh::global_n_processors() == 1)                                      \
        print_trace();                                                              \
      else                                                                          \
        libMesh::write_traceout();                                                  \
      libmesh_here();                                                               \
      MOOSE_ABORT;                                                                  \
    }                                                                               \
  } while (0)

/**
 * MooseWarning(), mooseInfo(), and mooseDeprecated() all print to _console instead of one of the
 * standard or wrapped streams directly. _console is defined as a true Global variable in Moose.h
 * so it's available anywhere these macros can be used. However, most of the time the macros are
 * called from a MooseObject-derived object where MooseObject::_console is available in a more
 * local scope. The one place where we've found these macros can break down is using them in a
 * static member function in a MooseObject derived class. In that case, the scope lookup finds
 * the MooseObject::_console object which is _not_ static instead of falling back to the global
 * version. This will cause a compile warning in that instance. The work-around for that rare
 * case is to print to the stream directly (not recommended), or make the method a const instance
 * method instead.
 */
#define mooseWarning(msg)                                                           \
  do                                                                                \
  {                                                                                 \
    if (Moose::_warnings_are_errors)                                                \
      mooseError(msg);                                                              \
    else                                                                            \
    {                                                                               \
      std::ostringstream _warn_oss_;                                                \
                                                                                    \
      _warn_oss_                                                                    \
        << COLOR_YELLOW                                                             \
        << "\n\n*** Warning ***\n"                                                  \
        << msg                                                                      \
        << "\nat " << __FILE__ << ", line " << __LINE__                             \
        << COLOR_DEFAULT                                                            \
        << "\n\n";                                                                  \
      if (Moose::_throw_on_error)                                                   \
        throw std::runtime_error(_warn_oss_.str());                                 \
      else                                                                          \
        _console << _warn_oss_.str() << std::flush;                                 \
    }                                                                               \
  } while (0)

#define mooseInfo(msg)                                                              \
    do                                                                              \
    {                                                                               \
      mooseDoOnce(                                                                  \
          _console                                                                  \
            << COLOR_CYAN                                                           \
            << "\n\n*** Info ***\n"                                                 \
            << msg                                                                  \
            << "\nat " << __FILE__ << ", line " << __LINE__                         \
            << COLOR_DEFAULT                                                        \
            << "\n" << std::endl;                                                   \
      );                                                                            \
    } while (0)

#define mooseDeprecated(msg)                                                                                  \
  do                                                                                                          \
  {                                                                                                           \
    if (Moose::_deprecated_is_error)                                                                          \
      mooseError("\n\nDeprecated code:\n" << msg << '\n');                                                    \
    else                                                                                                      \
      mooseDoOnce(                                                                                            \
        std::ostringstream _deprecate_oss_;                                                                   \
                                                                                                              \
        _deprecate_oss_                                                                                       \
          << COLOR_YELLOW                                                                                     \
          << "*** Warning, This code is deprecated, and likely to be removed in future library versions!\n"   \
          << msg << '\n'                                                                                      \
          << __FILE__ << ", line " << __LINE__ << ", compiled "                                               \
          << LIBMESH_DATE << " at " << LIBMESH_TIME << " ***"                                                 \
          << COLOR_DEFAULT                                                                                    \
          << "\n\n";                                                                                          \
        if (libMesh::global_n_processors() == 1)                                                              \
          print_trace(_deprecate_oss_);                                                                       \
        else                                                                                                  \
          libMesh::write_traceout();                                                                          \
        _console << _deprecate_oss_.str() << std::flush;                                                      \
      );                                                                                                      \
   } while (0)

#endif /* MOOSEERRORS_H */
