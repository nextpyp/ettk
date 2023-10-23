/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkDomainThreader_hxx
#define __itkDomainThreader_hxx

#include "itkDomainThreader.h"

namespace itk
{

template< class TDomainPartitioner, class TAssociate >
DomainThreader< TDomainPartitioner, TAssociate >
::DomainThreader()
{
  this->m_DomainPartitioner   = DomainPartitionerType::New();
  this->m_MultiThreader       = MultiThreader::New();
  this->m_NumberOfThreadsUsed = 0;
  this->m_Associate           = NULL;
}

template< class TDomainPartitioner, class TAssociate >
DomainThreader< TDomainPartitioner, TAssociate >
::~DomainThreader()
{
}

template< class TDomainPartitioner, class TAssociate >
MultiThreader *
DomainThreader< TDomainPartitioner, TAssociate >
::GetMultiThreader() const
{
  return this->m_MultiThreader;
}

template< class TDomainPartitioner, class TAssociate >
ThreadIdType
DomainThreader< TDomainPartitioner, TAssociate >
::GetMaximumNumberOfThreads() const
{
  return this->m_MultiThreader->GetNumberOfThreads();
}

template< class TDomainPartitioner, class TAssociate >
void
DomainThreader< TDomainPartitioner, TAssociate >
::SetMaximumNumberOfThreads( const ThreadIdType threads )
{
  if( threads != this->GetMaximumNumberOfThreads() )
    {
    this->m_MultiThreader->SetNumberOfThreads( threads );
    this->Modified();
    }
}

template< class TDomainPartitioner, class TAssociate >
void
DomainThreader< TDomainPartitioner, TAssociate >
::Execute( TAssociate * enclosingClass, const DomainType & completeDomain )
{
  this->m_Associate = enclosingClass;
  this->m_CompleteDomain = completeDomain;

  this->DetermineNumberOfThreadsUsed();

  this->BeforeThreadedExecution();

  // This calls ThreadedExecution in each thread.
  this->StartThreadingSequence();

  this->AfterThreadedExecution();
}

template< class TDomainPartitioner, class TAssociate >
void
DomainThreader< TDomainPartitioner, TAssociate >
::DetermineNumberOfThreadsUsed()
{
  DomainType subdomain;
  this->m_NumberOfThreadsUsed = this->m_DomainPartitioner->PartitionDomain(0,
                                            this->GetMultiThreader()->GetNumberOfThreads(),
                                            this->m_CompleteDomain,
                                            subdomain);
}

template< class TDomainPartitioner, class TAssociate >
void
DomainThreader< TDomainPartitioner, TAssociate >
::StartThreadingSequence()
{
  // Set up the multithreaded processing
  ThreadStruct str;
  str.domainThreader = this;

  MultiThreader* multiThreader = this->GetMultiThreader();
  multiThreader->SetSingleMethod(this->ThreaderCallback, &str);

  // multithread the execution
  multiThreader->SingleMethodExecute();
}

template< class TDomainPartitioner, class TAssociate >
ITK_THREAD_RETURN_TYPE
DomainThreader< TDomainPartitioner, TAssociate >
::ThreaderCallback( void* arg )
{
  MultiThreader::ThreadInfoStruct* info = static_cast<MultiThreader::ThreadInfoStruct *>(arg);
  ThreadStruct *str = static_cast<ThreadStruct *>(info->UserData);
  const ThreadIdType threadId    = info->ThreadID;
  const ThreadIdType threadCount = info->NumberOfThreads;

  // Get the sub-domain to process for this thread.
  DomainType subdomain;
  const ThreadIdType total = str->domainThreader->GetDomainPartitioner()->PartitionDomain(threadId,
                                            threadCount,
                                            str->domainThreader->m_CompleteDomain,
                                            subdomain);

  // Execute the actual method with appropriate sub-domain.
  // If the threadID is greater than the total number of regions
  // that PartitionDomain will create, don't use this thread.
  // Sometimes the threads dont break up very well and it is just
  // as efficient to leave a few threads idle.
  if ( threadId < total )
    {
    str->domainThreader->ThreadedExecution( subdomain, threadId );
    }

  return ITK_THREAD_RETURN_VALUE;
}
}

#endif
