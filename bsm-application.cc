/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2014 North Carolina State University
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Scott E. Carpenter <scarpen@ncsu.edu>
 *
 */

#include "../../../build/ns3/bsm-application.h"
#include "bsm-application.h"


#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "../../../build/ns3/address.h"
#include "../../../build/ns3/assert.h"
#include "../../../build/ns3/buffer.h"
#include "../../../build/ns3/callback.h"
#include "../../../build/ns3/inet-socket-address.h"
#include "../../../build/ns3/ipv4.h"
#include "../../../build/ns3/ipv4-address.h"
#include "../../../build/ns3/ipv4-interface-container.h"
#include "../../../build/ns3/log.h"
#include "../../../build/ns3/log-macros-disabled.h"
#include "../../../build/ns3/mobility-helper.h"
#include "../../../build/ns3/mobility-model.h"
#include "../../../build/ns3/net-device.h"
#include "../../../build/ns3/node.h"
#include "../../../build/ns3/nstime.h"
#include "../../../build/ns3/object.h"
#include "../../../build/ns3/object-base.h"
#include "../../../build/ns3/packet.h"
#include "../../../build/ns3/ptr.h"
#include "../../../build/ns3/simulator.h"
#include "../../../build/ns3/socket.h"
#include "../../../build/ns3/type-id.h"
#include "ns3/qos-tag.h"

NS_LOG_COMPONENT_DEFINE ("BsmApplicationMing");

namespace ns3 {

// (Arbitrary) port for establishing socket to transmit WAVE BSMs
int BsmApplication::wavePort = 9080;

NS_OBJECT_ENSURE_REGISTERED (BsmApplication);


/*Ming*/
FwdHeader::FwdHeader (uint16_t fwd)
  : m_fwd (fwd)
{
}

FwdHeader::~FwdHeader ()
{
}

TypeId
FwdHeader::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::FwdHeader")
    .SetParent<Header> ()
    .SetGroupName ("Fwd")
    .AddConstructor<FwdHeader> ();
  return tid;
}

TypeId
FwdHeader::GetInstanceTypeId () const
{
  return GetTypeId ();
}

uint32_t
FwdHeader::GetSerializedSize () const
{
  return 2;
}

void
FwdHeader::Serialize (Buffer::Iterator i) const
{
  i.WriteHtonU16 (m_fwd);
}

uint32_t
FwdHeader::Deserialize (Buffer::Iterator start)
{
  Buffer::Iterator i = start;

  m_fwd = i.ReadNtohU16 ();

  uint32_t dist = i.GetDistanceFrom (start);
  NS_ASSERT (dist == GetSerializedSize ());
  return dist;
}

void
FwdHeader::Print (std::ostream &os) const
{
  os << "HasItBeForwarded: " << m_fwd;
}
/*Ming*/


TypeId
BsmApplication::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::BsmApplication")
    .SetParent<Application> ()
    .SetGroupName ("Wave")
    .AddConstructor<BsmApplication> ()
    ;
  return tid;
}

BsmApplication::BsmApplication ()
  : m_waveBsmStats (0),
    m_txSafetyRangesSq (),
    m_TotalSimTime (Seconds (10)),
    m_wavePacketSize (200),
    m_numWavePackets (1),
    m_waveInterval (MilliSeconds (100)),
    ForwarderNode (false),
    m_gpsAccuracyNs (10000),
    m_adhocTxInterfaces (0),
    m_nodesMoving (0),
    m_unirv (0),
    m_nodeId (0),
    m_chAccessMode (0),
    m_txMaxDelay (MilliSeconds (10)), //10 msec
    m_prevTxDelay (MilliSeconds (0))
{
  NS_LOG_FUNCTION (this);
}

BsmApplication::~BsmApplication ()
{
  NS_LOG_FUNCTION (this);
}

void
BsmApplication::DoDispose (void)
{
  NS_LOG_FUNCTION (this);

  // chain up
  Application::DoDispose ();
}

// Application Methods
void BsmApplication::StartApplication () // Called at time specified by Start
{
  NS_LOG_FUNCTION (this);


/*Ming*/
 // if((m_node->GetId()==33)||(m_node->GetId()==35)||(m_node->GetId()==9)||(m_node->GetId()==3)||(m_node->GetId()==15)||(m_node->GetId()==29)||(m_node->GetId()==3)||(m_node->GetId()==15)||(m_node->GetId()==29)||(m_node->GetId()==22)||(m_node->GetId()==32)||(m_node->GetId()==17))
  	 // if((m_node->GetId()==3))
  //  if((m_node->GetId()==33)||(m_node->GetId()==35)||(m_node->GetId()==9))
	//if((m_node->GetId()==33)||(m_node->GetId()==35)||(m_node->GetId()==9)||(m_node->GetId()==3)||(m_node->GetId()==15)||(m_node->GetId()==29))
	//if((m_node->GetId()==33)||(m_node->GetId()==35)||(m_node->GetId()==9)||(m_node->GetId()==3)||(m_node->GetId()==15)||(m_node->GetId()==29)||(m_node->GetId()==12)||(m_node->GetId()==15)||(m_node->GetId()==8))
  	  //if((m_node->GetId()==33)||(m_node->GetId()==35)||(m_node->GetId()==9)||(m_node->GetId()==3)||(m_node->GetId()==15)||(m_node->GetId()==29)||(m_node->GetId()==12)||(m_node->GetId()==15)||(m_node->GetId()==8)||(m_node->GetId()==16)||(m_node->GetId()==4)||(m_node->GetId()==39))
	 // if((m_node->GetId()==33)||(m_node->GetId()==35)||(m_node->GetId()==9)||(m_node->GetId()==3)||(m_node->GetId()==15)||(m_node->GetId()==29)||(m_node->GetId()==12)||(m_node->GetId()==15)||(m_node->GetId()==8)||(m_node->GetId()==16)||(m_node->GetId()==4)||(m_node->GetId()==39)||(m_node->GetId()==5)||(m_node->GetId()==7)||(m_node->GetId()==17))
  {
	  //SetForwarderNode(false);
	 SetForwarderNode(true);
	  NS_LOG_UNCOND("Node "<<m_node->GetId());
	  //  SetForwarderNode(true);
  }
/*Ming*/

	 // setup generation of WAVE BSM messages
  Time waveInterPacketInterval = m_waveInterval;

  // BSMs are not transmitted for the first second
  Time startTime = Seconds (1.0); //1.0 default
  // total length of time transmitting WAVE packets
  Time totalTxTime = m_TotalSimTime - startTime;
  // total WAVE packets needing to be sent
  m_numWavePackets = (uint32_t) (totalTxTime.GetDouble () / m_waveInterval.GetDouble ());

  TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");

  // every node broadcasts WAVE BSM to potentially all other nodes
  Ptr<Socket> recvSink = Socket::CreateSocket (GetNode (m_nodeId), tid);
  recvSink->SetRecvCallback (MakeCallback (&BsmApplication::ReceiveWavePacket, this));
  InetSocketAddress local = InetSocketAddress (Ipv4Address::GetAny (), wavePort);
  recvSink->Bind (local);
  recvSink->BindToNetDevice (GetNetDevice (m_nodeId));
  recvSink->SetAllowBroadcast (true);

  // dest is broadcast address
  InetSocketAddress remote = InetSocketAddress (Ipv4Address ("255.255.255.255"), wavePort);
  recvSink->Connect (remote);

  // Transmission start time for each BSM:
  // We assume that the start transmission time
  // for the first packet will be on a ns-3 time
  // "Second" boundary - e.g., 1.0 s.
  // However, the actual transmit time must reflect
  // additional effects of 1) clock drift and
  // 2) transmit delay requirements.
  // 1) Clock drift - clocks are not perfectly
  // synchronized across all nodes.  In a VANET
  // we assume all nodes sync to GPS time, which
  // itself is assumed  accurate to, say, 40-100 ns.
  // Thus, the start transmission time must be adjusted
  // by some value, t_drift.
  // 2) Transmit delay requirements - The US
  // minimum performance requirements for V2V
  // BSM transmission expect a random delay of
  // +/- 5 ms, to avoid simultanous transmissions
  // by all vehicles congesting the channel.  Thus,
  // we need to adjust the start trasmission time by
  // some value, t_tx_delay.
  // Therefore, the actual transmit time should be:
  // t_start = t_time + t_drift + t_tx_delay
  // t_drift is always added to t_time.
  // t_tx_delay is supposed to be +/- 5ms, but if we
  // allow negative numbers the time could drift to a value
  // BEFORE the interval start time (i.e., at 100 ms
  // boundaries, we do not want to drift into the
  // previous interval, such as at 95 ms.  Instead,
  // we always want to be at the 100 ms interval boundary,
  // plus [0..10] ms tx delay.
  // Thus, the average t_tx_delay will be
  // within the desired range of [0..10] ms of
  // (t_time + t_drift)

  // WAVE devices sync to GPS time
  // and all devices would like to begin broadcasting
  // their safety messages immediately at the start of
  // the CCH interval.  However, if all do so, then
  // significant collisions occur.  Thus, we assume there
  // is some GPS sync accuracy on GPS devices,
  // typically 40-100 ns.
  // Get a uniformly random number for GPS sync accuracy, in ns.
  Time tDrift = NanoSeconds (m_unirv->GetInteger (0, m_gpsAccuracyNs));

  // When transmitting at a default rate of 10 Hz,
  // the subsystem shall transmit every 100 ms +/-
  // a random value between 0 and 5 ms. [MPR-BSMTX-TXTIM-002]
  // Source: CAMP Vehicle Safety Communications 4 Consortium
  // On-board Minimum Performance Requirements
  // for V2V Safety Systems Version 1.0, December 17, 2014
  // max transmit delay (default 10ms)
  // get value for transmit delay, as number of ns
  uint32_t d_ns = static_cast<uint32_t> (m_txMaxDelay.GetInteger ());
  // convert random tx delay to ns-3 time
  // see note above regarding centering tx delay
  // offset by 5ms + a random value.
  Time txDelay = NanoSeconds (m_unirv->GetInteger (0, d_ns));
  m_prevTxDelay = txDelay;

  Time txTime = startTime + tDrift + txDelay;
 NS_LOG_UNCOND("txTime= "<<txTime);
 NS_LOG_UNCOND("Simu Time now="<<Simulator::Now().GetNanoSeconds()<<" startTime= "<<startTime);
 // schedule transmission of first packet
  Simulator::ScheduleWithContext (recvSink->GetNode ()->GetId (),
                                  txTime, &BsmApplication::GenerateWaveTraffic, this,
                                  recvSink, m_wavePacketSize, m_numWavePackets, waveInterPacketInterval, m_nodeId);
}

void BsmApplication::StopApplication () // Called at time specified by Stop
{
  NS_LOG_FUNCTION (this);
}

/*Ming*/
void
BsmApplication::SetForwarderNode(bool f)
{
	ForwarderNode=f;
}

/*bool
BsmApplication::GetForwarderNode()
{
	return ForwarderNode;
}*/
/*Ming*/


bool
BsmApplication::GetForwarderNode(Ptr<Node> rxNodeCurr)
{

	NS_LOG_FUNCTION (this);
	//Initialization
	   double MAX= 0.0;
	   int rxNodeTargetId;

       // Find which candidate node is furthest
	   // After finding rxNodeTargetId
	   // compare between rxNodeCurr and rxNodeTargetId, to decide whether
	   //rxNodeCurr is the forwarding not or not
       int nRxNodes = m_adhocTxInterfaces->GetN ();
       int nTxNodes = m_adhocTxInterfaces->GetN ();
       for (int i = 0; i < nRxNodes; i++)
         {
           Ptr<Node> rxNode = GetNode (i); //receiver node
           int rxNodeId = rxNode->GetId (); //receiver node id
           for(int j=0; j<nTxNodes; j++) //loop for comparing with all transmitting nodes
            {
           		 Ptr<Node> txNode = GetNode (j);
           		 int txNodeId = txNode->GetId();


           		 if (rxNodeId != txNodeId)
           		 {
           			 Ptr<MobilityModel> rxPosition = rxNode->GetObject<MobilityModel> ();
           			 NS_ASSERT_MSG (rxPosition != 0, "Position error");
               // confirm that the receiving node
               // has also started moving in the scenario
               // if it has not started moving, then
               // it is not a candidate to receive a packet
           			 int receiverMoving = m_nodesMoving->at (rxNodeId);
           			 if (receiverMoving == 1)
           			 {
           				 double distSq = MobilityHelper::GetDistanceSquaredBetween (txNode, rxNode);
           				 double dist = sqrt(distSq); //calculate the distance
           				//double distance = rxPosition->GetDistanceFrom (txPosition); //Also possible to calculate distance in this way
           				 // dest node within range?
           				 if (distSq > 0.0 && dist<=300) //checking only for 150 transmission range of  a receiver
           				 {
           					 if(MAX<dist)
           					 {
           						 MAX=dist;
           						 rxNodeTargetId=rxNodeId;
           					 	 }//end for MAX

           				 	 }//end if (distSq>0.0)
           			 	 }//end  if (receiverMoving == 1)
           		 	 }//end  if (rxNodeId != txNodeId)

           	 	 }//end inner for loop
         }//Outer for loop


      // std::cout<<"RcverNodeID"<<rxNodeCurr->GetId()<<" TargetNodeID"<<rxNodeTargetId<<" Max Distance"<<MAX;

      int rxNodeCurrId = rxNodeCurr->GetId();
       if(rxNodeCurrId == rxNodeTargetId)
           	   return true;
    	   else
    		   return false;

}//end GetForwarderNode(Ptr<Node> rxNodeCurr)

/*End Nawaz*/




void
BsmApplication::Setup (Ipv4InterfaceContainer & i,
                       int nodeId,
                       Time totalTime,
                       uint32_t wavePacketSize, // bytes
                       Time waveInterval,
                       double gpsAccuracyNs,
                       std::vector <double> rangesSq,           // m ^2
                       Ptr<WaveBsmStats> waveBsmStats,
                       std::vector<int> * nodesMoving,
                       int chAccessMode,
                       Time txMaxDelay)
{
  NS_LOG_FUNCTION (this);

  m_unirv = CreateObject<UniformRandomVariable> ();

  m_TotalSimTime = totalTime;
  m_wavePacketSize = wavePacketSize;
  m_waveInterval = waveInterval;
  m_gpsAccuracyNs = gpsAccuracyNs;
  int size = rangesSq.size ();
  m_waveBsmStats = waveBsmStats;
  m_nodesMoving = nodesMoving;
  m_chAccessMode = chAccessMode;
  m_txSafetyRangesSq.clear ();
  m_txSafetyRangesSq.resize (size, 0);

  for (int index = 0; index < size; index++)
    {
      // stored as square of value, for optimization
      m_txSafetyRangesSq[index] = rangesSq[index];
    }

  m_adhocTxInterfaces = &i;
  m_nodeId = nodeId;
  m_txMaxDelay = txMaxDelay;
}

/*Ming*/
void
BsmApplication::PrintNodeInfo ()
{
	NS_LOG_FUNCTION (this);
//	if(m_nodeId==0)
	//NS_LOG_UNCOND ("Node "<< m_nodeId << "Number of Packets transmitted "<< m_waveBsmStats->GetTxPktCount());
	//NS_LOG_UNCOND ("Node "<< m_nodeId << "Number of Packets received "<< m_waveBsmStats->GetRxPktCount());

}
/*Ming*/

void
BsmApplication::GenerateWaveTraffic (Ptr<Socket> socket, uint32_t pktSize,
                                     uint32_t pktCount, Time pktInterval,
                                     uint32_t sendingNodeId)
{
	NS_LOG_FUNCTION (this);

	//NS_LOG_FUNCTION ("Node "<<m_node->GetId()<<" GenerateWaveTraffic");

	//NS_LOG_UNCOND("Time now="<<Simulator::Now().GetNanoSeconds());
	//NS_LOG_UNCOND("pktCount= "<<pktCount );
	// more packets to send?

	//Time nodeDelay = NanoSeconds (m_unirv->GetInteger (sendingNodeId, 30));

	if (pktCount > 0)
    {
      // for now, we cannot tell if each node has
      // started mobility.  so, as an optimization
      // only send if  this node is moving
      // if not, then skip
      int txNodeId = sendingNodeId;
      Ptr<Node> txNode = GetNode (txNodeId);
      Ptr<MobilityModel> txPosition = txNode->GetObject<MobilityModel> ();
      Vector pos_Tx= txPosition->GetPosition();

      //NS_LOG_UNCOND("txnodeID="<<txNodeId<<"pktCount= "<<pktCount);
     // if(txNodeId==0)
    //	  NS_LOG_UNCOND("txnodeID="<<txNodeId<<"time now="<<Simulator::Now().GetMicroSeconds());
      //NS_ASSERT (txPosition != 0); //Without mobility checking nawaz

      int senderMoving = m_nodesMoving->at (txNodeId);
     // if (senderMoving != 0) //Without mobility checking nawaz
     // {
    	  Ptr<Packet> p=Create<Packet> (pktSize);

    	  /*Ming*/
		  FwdHeader fwdheader;
		  fwdheader.SetFwd(0);
    	  p->AddHeader(fwdheader) ;
    	  	//* End Ming *//


    	  /*AC_BE = 0; AC_BK = 1; AC_VI = 2;AC_VO = 3*/

    	  uint8_t qos_tid;
    	  if(pktCount%4 == 0)
    		   qos_tid = 0;
    	  else if(pktCount%4 == 1)
    		   qos_tid = 1;
    	  else if(pktCount%4 == 2)
    	       qos_tid = 2;
    	  else if(pktCount%4 == 3)
    	       qos_tid = 3;

    	  QosTag qosTag;
    	  qosTag.SetTid(qos_tid);
    	  p->AddPacketTag(qosTag);

    	  // send it!
 //         socket->Send (Create<Packet> (pktSize));
          socket->Send (p);
          /*Ming*/
          // count it
         // NS_LOG_UNCOND("pkt size="<<p->GetSize());
          /*Transmitter position*/
          //if(pos_Tx.x>=-10 && pos_Tx.x<=10) //transmitter condition for position 'G'
          // if(pos_Tx.x>=-40 && pos_Tx.x<=-10) //transmitter condition for position 'A1'
           if(pos_Tx.x<-70 && pos_Tx.x>=-150) //transmitter condition for position 'A3'
          {
          m_waveBsmStats->IncTxPktCount ();
          m_waveBsmStats->IncTxByteCount (pktSize);
          int wavePktsSent = m_waveBsmStats->GetTxPktCount ();
          if ((m_waveBsmStats->GetLogging () != 0) && ((wavePktsSent % 1000) == 0))
            {
              NS_LOG_UNCOND ("Sending WAVE pkt # " << wavePktsSent );
            }

          // find other nodes within range that would be
          // expected to receive this broadcast
          int nRxNodes = m_adhocTxInterfaces->GetN ();
          for (int i = 0; i < nRxNodes; i++)
            {
              Ptr<Node> rxNode = GetNode (i);
              int rxNodeId = rxNode->GetId ();

              if (rxNodeId != txNodeId)
                {
                  Ptr<MobilityModel> rxPosition = rxNode->GetObject<MobilityModel> ();
                  Vector pos_Rx = rxPosition->GetPosition();
                  //NS_ASSERT (rxPosition != 0); //Without mobility checking nawaz
                  // confirm that the receiving node
                  // has also started moving in the scenario
                  // if it has not started moving, then
                  // it is not a candidate to receive a packet
                  int receiverMoving = m_nodesMoving->at (rxNodeId);
                 // if (receiverMoving == 1)
                /*Receiver position*/
                   if(pos_Rx.x>=-10 && pos_Rx.x<=10)//receiver at position 'G'
                //if((pos_Rx.x>=-70 && pos_Rx.x<-10)||(pos_Rx.x>10 && pos_Rx.x<=70))//receiver at position 'A1/A2/B1/B2'
               // if(pos_Rx.x>=-70 && pos_Rx.x<-10)//receiver at position 'A1/A2
                //if((pos_Rx.x<-70 && pos_Rx.x>=-150)||(pos_Rx.x>70 && pos_Rx.x<=150))//receiver at position 'A3/B3'
              //  if(pos_Rx.x<-70 && pos_Rx.x>=-150)//receiver at position 'A3'
                // if(pos_Rx.x>10 && pos_Rx.x<=40)//receiver at position 'B1'
                 {
                      double distSq = MobilityHelper::GetDistanceSquaredBetween (txNode, rxNode);
                      if (distSq > 0.0)
                        {
                          // dest node within range?
                          int rangeCount = m_txSafetyRangesSq.size ();
                          for (int index = 1; index <= rangeCount; index++)
                            {
                              if (distSq <= m_txSafetyRangesSq[index - 1])
                                {
                                  // we should expect dest node to receive broadcast pkt
                                  m_waveBsmStats->IncExpectedRxPktCount (index);
                                }
                            }
                        }
                    }
                }
            }
        }

      // every BSM must be scheduled with a tx time delay
      // of +/- (5) ms.  See comments in StartApplication().
      // we handle this as a tx delay of [0..10] ms
      // from the start of the pktInterval boundary
      uint32_t d_ns = static_cast<uint32_t> (m_txMaxDelay.GetInteger ());
      Time txDelay = NanoSeconds (m_unirv->GetInteger (0, d_ns));
     // Time tDrift = NanoSeconds (m_unirv->GetInteger (0, m_gpsAccuracyNs)); //added by Nawaz
     // Time tDrift = NanoSeconds (m_unirv->GetValue(0.0, m_gpsAccuracyNs));
      // do not want the tx delay to be cumulative, so
      // deduct the previous delay value.  thus we adjust
      // to schedule the next event at the next pktInterval,
      // plus some new [0..10] ms tx delay
      //Time txTime = pktInterval - m_prevTxDelay + txDelay; // + nodeDelay;
      //m_prevTxDelay = txDelay;

      Time txTime = pktInterval + txDelay;//+tDrift; //tDrift added by Nawaz

      /**Nawaz**/
      if(txTime<0.0)
      {
    	  txTime=txDelay;
    	  m_prevTxDelay=txDelay;
      }
    	  /**Nawaz**/

      //NS_LOG_UNCOND("Time Now= "<<Simulator::Now().GetNanoSeconds()<<" pktInterval= "<<pktInterval<<" txTime = "<<txTime<<"txDelay= "<<txDelay);

      Simulator::ScheduleWithContext (socket->GetNode ()->GetId (),
                                      txTime, &BsmApplication::GenerateWaveTraffic, this,
                                      socket, pktSize, pktCount - 1, pktInterval,  socket->GetNode ()->GetId ());
    }
  else
    {
      socket->Close ();
    }
}

void BsmApplication::ReceiveWavePacket (Ptr<Socket> socket)
{
	//NS_LOG_FUNCTION ("Node "<<m_node->GetId()<<" ReceiveWavePacket");//Ming

  Ptr<Packet> packet;
  while ((packet = socket->Recv ()))
    {
      Ptr<Node> rxNode = socket->GetNode ();

  /*Ming*/
      // 	  std::ostringstream oss;    packet->Print (oss);    NS_LOG_UNCOND( oss.str ());

   	 // NS_LOG_UNCOND("rcv pkt size="<<packet->GetSize());
      FwdHeader fwdheader;
      packet->RemoveHeader(fwdheader);
      NS_LOG_FUNCTION ("Packet ID "<<packet->GetUid());
      /*Ming*/

      SocketAddressTag tag;
      bool found;
      found = packet->PeekPacketTag (tag);

      if (found)
        {

          InetSocketAddress addr = InetSocketAddress::ConvertFrom (tag.GetAddress ());
          int nodes = m_adhocTxInterfaces->GetN ();
          for (int i = 0; i < nodes; i++)
            {
        	  if (addr.GetIpv4 () == m_adhocTxInterfaces->GetAddress (i) )
        	  {
              Ptr<Node> txNode = GetNode (i);

              /*Ming rebroadcasting is done here*/
             /* //Apply your logic in GetForwarderNode() function, the node is the furthest will rebroadcast, otherwise will not do rebroadcast
             if(fwdheader.GetFwd() == 0 && GetForwarderNode(rxNode))//Nawaz
            //if(fwdheader.GetFwd() == 0 && rxNode->GetId()== 23)//Nawaz rxNode->GetId()== 3, only one particular node as a rebroadcaster such as nodeID 3
            //	if(fwdheader.GetFwd() == 0 && GetForwarderNode())//Ming
              {
            	 NS_LOG_UNCOND(" ForwarderNodeID = "<<rxNode->GetId());
            	 packet->RemovePacketTag(tag);
            	  ///Invoke for relaying mechanism here//
            	  RelayReceivedBsmPacket(socket, packet, packet->GetSize(), rxNode->GetId()); //Rebrodcast the BSM
              } */
             /*Ming*/

             	 	 //HandleReceivedInvoBSMPacket(packet, txNode, rxNode);//Ming (actually done thing, just showing some statistics)
            //  if(txNode->GetId()==0)
              //    NS_LOG_UNCOND("Rcv txnodeID="<<txNode->GetId()<<"time now="<<Simulator::Now().GetMicroSeconds());

              Ptr<MobilityModel> txPosition = txNode->GetObject<MobilityModel> ();
                    Vector pos_Tx= txPosition->GetPosition();

                    Ptr<MobilityModel> rxPosition = rxNode->GetObject<MobilityModel>();
                    Vector pos_Rx = rxPosition->GetPosition();

                   //Transmitter position condition
                   //if(pos_Tx.x>=-10 && pos_Tx.x<=10) //transmitter condition for position 'G'
                   // if(pos_Tx.x>=-40 && pos_Tx.x<-10) //transmitter condition for position 'A1'
                    if(pos_Tx.x<-70 && pos_Tx.x>=-150) //transmitter condition for position 'A3'
                    /*Receiver position condition*/
                    {
                    	/*Receiver position condition*/
                    if(pos_Rx.x>=-10 && pos_Rx.x<=10) //receiver at position 'G'
                  //  if((pos_Rx.x>=-70 && pos_Rx.x<-10)||(pos_Rx.x>10 && pos_Rx.x<=70))//receiver at position 'A1/A2/B1/B2'
                   // if(pos_Rx.x>=-70 && pos_Rx.x<-10) //receiver at position A1/A2
                    	// if((pos_Rx.x<-70 && pos_Rx.x>=-150)||(pos_Rx.x>70 && pos_Rx.x<=150))//receiver at position 'A3/B3'
                    // if(pos_Rx.x<-70 && pos_Rx.x>=-150)//receiver at position 'A3'
                   //  if(pos_Rx.x>10 && pos_Rx.x<=40)//receiver at position 'B1'
                    	HandleReceivedBsmPacket (txNode, rxNode);//collection BSMs receive stats
                    }

        	  }
            }
        }
    }
}

/*RelayReceivedBsmPacket added by Ming*/
void
BsmApplication::RelayReceivedBsmPacket (Ptr<Socket> socket, Ptr<Packet> p, uint32_t pktSize, uint32_t sendingNodeId)
{
	//NS_LOG_FUNCTION ("Node "<<m_node->GetId()<<" RelayReceivedBsmPacket; "<<"Packet ID "<<p->GetUid());

	int txNodeId = sendingNodeId;
	Ptr<Node> txNode = GetNode (txNodeId);
	Ptr<MobilityModel> txPosition = txNode->GetObject<MobilityModel> ();
	NS_ASSERT (txPosition != 0);

	int senderMoving = m_nodesMoving->at (txNodeId);

	if (senderMoving != 0)
	{
		FwdHeader fwdheader;
		fwdheader.SetFwd(1);
  	    p->AddHeader(fwdheader) ;
		socket->Send (p);

		m_waveBsmStats->IncTxPktCount ();
		          m_waveBsmStats->IncTxByteCount (pktSize);
		          int wavePktsSent = m_waveBsmStats->GetTxPktCount ();
		          if ((m_waveBsmStats->GetLogging () != 0) && ((wavePktsSent % 1000) == 0))
		            {
		              NS_LOG_UNCOND ("Sending WAVE pkt # " << wavePktsSent );
		            }

		          // find other nodes within range that would be
		          // expected to receive this broadbast
		          int nRxNodes = m_adhocTxInterfaces->GetN ();
		          for (int i = 0; i < nRxNodes; i++)
		            {
		              Ptr<Node> rxNode = GetNode (i);
		              int rxNodeId = rxNode->GetId ();

		              if (rxNodeId != txNodeId)
		                {
		                  Ptr<MobilityModel> rxPosition = rxNode->GetObject<MobilityModel> ();
		                  NS_ASSERT (rxPosition != 0);
		                  // confirm that the receiving node
		                  // has also started moving in the scenario
		                  // if it has not started moving, then
		                  // it is not a candidate to receive a packet
		                  int receiverMoving = m_nodesMoving->at (rxNodeId);
		                  if (receiverMoving == 1)
		                    {
		                      double distSq = MobilityHelper::GetDistanceSquaredBetween (txNode, rxNode);
		                      if (distSq > 0.0)
		                        {
		                          // dest node within range?
		                          int rangeCount = m_txSafetyRangesSq.size ();
		                          for (int index = 1; index <= rangeCount; index++)
		                            {
		                              if (distSq <= m_txSafetyRangesSq[index - 1])
		                                {
		                                  // we should expect dest node to receive broadcast pkt
		                                  m_waveBsmStats->IncExpectedRxPktCount (index);
		                                }
		                            }
		                        }
		                    }
		                }
		            }
		        }
}

/*ReceiveRelayedBsm added by Ming*/
void
BsmApplication::ReceiveRelayedBsm(Ptr<Socket> socket)
{

}
/*HandleReceivedInvoBSMPacket added by Ming*/
void
BsmApplication::HandleReceivedInvoBSMPacket (Ptr<Packet> p,Ptr<Node> txNode, Ptr<Node> rxNode )
{
	NS_LOG_FUNCTION (this);
	  Ptr<MobilityModel> rxPosition = rxNode->GetObject<MobilityModel> ();
	  NS_ASSERT (rxPosition != 0);
	  // confirm that the receiving node
	  // has also started moving in the scenario
	  // if it has not started moving, then
	  // it is not a candidate to receive a packet
	  int rxNodeId = rxNode->GetId ();
	  int receiverMoving = m_nodesMoving->at (rxNodeId);
	  if (receiverMoving == 1)
	    {
	      double rxDistSq = MobilityHelper::GetDistanceSquaredBetween (rxNode, txNode);
	      if (rxDistSq > 0.0)
	        {
	          int rangeCount = m_txSafetyRangesSq.size ();
	          for (int index = 1; index <= rangeCount; index++)
	            {
	              if (rxDistSq <= m_txSafetyRangesSq[index - 1])
	              {

	            	 // NS_LOG_UNCOND("Node "<< m_nodeId<<" pktID "<<p->GetUid());
	            	 // std::ostringstream oss;    p->Print (oss);    NS_LOG_UNCOND( oss.str ());
	            	  //if(m_waveBsmStats->RxInocheck(index, p->GetUid()))
	                  //m_waveBsmStats->IncRxInvoInRangeCount (index);
	                }

	           //   m_waveBsmStats->PrintPktIDlist();
	            }
	        }
	    }

}

void BsmApplication::HandleReceivedBsmPacket (Ptr<Node> txNode,
                                              Ptr<Node> rxNode)
{
  NS_LOG_FUNCTION (this);

  m_waveBsmStats->IncRxPktCount ();

  Ptr<MobilityModel> rxPosition = rxNode->GetObject<MobilityModel> ();

  //NS_ASSERT (rxPosition != 0);//Without mobility checking nawaz
  // confirm that the receiving node
  // has also started moving in the scenario
  // if it has not started moving, then
  // it is not a candidate to receive a packet
  int rxNodeId = rxNode->GetId ();
  int receiverMoving = m_nodesMoving->at (rxNodeId);
  //if (receiverMoving == 1) //Without mobility checking Nawaz
    {
      double rxDistSq = MobilityHelper::GetDistanceSquaredBetween (rxNode, txNode);
      if (rxDistSq > 0.0)
        {
          int rangeCount = m_txSafetyRangesSq.size ();
          for (int index = 1; index <= rangeCount; index++)
            {
              if (rxDistSq <= m_txSafetyRangesSq[index - 1])
                {
                  m_waveBsmStats->IncRxPktInRangeCount (index);
                }
            }
        }
    }
}

int64_t
BsmApplication::AssignStreams (int64_t streamIndex)
{
  NS_LOG_FUNCTION (this);

  NS_ASSERT (m_unirv);  // should be set by Setup() prevoiusly
  m_unirv->SetStream (streamIndex);

  return 1;
}

Ptr<Node>
BsmApplication::GetNode (int id)
{
  NS_LOG_FUNCTION (this);

  std::pair<Ptr<Ipv4>, uint32_t> interface = m_adhocTxInterfaces->Get (id);
  Ptr<Ipv4> pp = interface.first;
  Ptr<Node> node = pp->GetObject<Node> ();

  return node;
}

Ptr<NetDevice>
BsmApplication::GetNetDevice (int id)
{
  NS_LOG_FUNCTION (this);

  std::pair<Ptr<Ipv4>, uint32_t> interface = m_adhocTxInterfaces->Get (id);
  Ptr<Ipv4> pp = interface.first;
  Ptr<NetDevice> device = pp->GetObject<NetDevice> ();

  return device;
}

} // namespace ns3
