-- MySQL dump 10.19  Distrib 10.3.32-MariaDB, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: ProjectLipoate2
-- ------------------------------------------------------
-- Server version	10.3.32-MariaDB-0ubuntu0.20.04.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8mb4 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `Clusters`
--

DROP TABLE IF EXISTS `Clusters`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Clusters` (
  `cluster_ID` varchar(32) NOT NULL,
  `genome_ID` varchar(32) DEFAULT NULL,
  PRIMARY KEY (`cluster_ID`),
  KEY `genome_ID` (`genome_ID`),
  CONSTRAINT `Clusters_ibfk_1` FOREIGN KEY (`genome_ID`) REFERENCES `Genomes` (`genome_ID`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Functions`
--

DROP TABLE IF EXISTS `Functions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Functions` (
  `ID` int(11) NOT NULL AUTO_INCREMENT,
  `cluster_ID` varchar(32) DEFAULT NULL,
  `keyword` varchar(32) DEFAULT NULL,
  PRIMARY KEY (`ID`),
  KEY `cluster_ID` (`cluster_ID`),
  CONSTRAINT `Functions_ibfk_1` FOREIGN KEY (`cluster_ID`) REFERENCES `Clusters` (`cluster_ID`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=159598 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Genomes`
--

DROP TABLE IF EXISTS `Genomes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Genomes` (
  `genome_ID` varchar(32) NOT NULL,
  `status` varchar(3) DEFAULT NULL,
  `Phylum` varchar(128) DEFAULT NULL,
  `Class` varchar(128) DEFAULT NULL,
  `Ordnung` varchar(128) DEFAULT NULL,
  `Family` varchar(128) DEFAULT NULL,
  `Genus` varchar(128) DEFAULT NULL,
  `Species` varchar(128) DEFAULT NULL,
  `strain` varchar(128) DEFAULT NULL,
  `TypeStrain` tinyint(4) DEFAULT NULL,
  `Completeness` decimal(5,2) DEFAULT NULL,
  `Contamination` decimal(5,2) DEFAULT NULL,
  `dRep` tinyint(1) DEFAULT NULL,
  `NCBITaxon` int(11) DEFAULT NULL,
  `NCBIProject` int(11) DEFAULT NULL,
  `NCBIBioproject` varchar(32) DEFAULT NULL,
  `NCBIBiosample` varchar(32) DEFAULT NULL,
  `NCBIAssembly` varchar(32) DEFAULT NULL,
  `Updated` datetime DEFAULT NULL,
  PRIMARY KEY (`genome_ID`),
  UNIQUE KEY `NCBIAssembly` (`NCBIAssembly`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Proteins`
--

DROP TABLE IF EXISTS `Proteins`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Proteins` (
  `protein_ID` varchar(32) NOT NULL,
  `genome_ID` varchar(32) DEFAULT NULL,
  `locustag` varchar(32) DEFAULT NULL,
  `contig` varchar(32) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  `end` int(11) DEFAULT NULL,
  `strand` varchar(1) DEFAULT NULL,
  `cluster_ID` varchar(32) DEFAULT NULL,
  `BBH` varchar(32) DEFAULT NULL,
  `PHY` varchar(32) DEFAULT NULL,
  `HMM` varchar(32) DEFAULT NULL,
  `score` smallint(6) DEFAULT NULL,
  `domains` smallint(6) DEFAULT NULL,
  `sequence` varchar(2048) DEFAULT NULL,
  PRIMARY KEY (`protein_ID`),
  KEY `genome_ID` (`genome_ID`),
  KEY `cluster_ID` (`cluster_ID`),
  CONSTRAINT `Proteins_ibfk_1` FOREIGN KEY (`genome_ID`) REFERENCES `Genomes` (`genome_ID`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `Proteins_ibfk_2` FOREIGN KEY (`cluster_ID`) REFERENCES `Clusters` (`cluster_ID`) ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2022-01-26 22:24:09
