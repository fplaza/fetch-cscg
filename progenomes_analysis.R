cscg.dir='/mnt/nfsdata/dev/florian/fetch-cscg/cscg_db'
cscg.bacteria=read.table(paste(cscg.dir,'cscg_bacteria.tab', sep='/'),header=T, stringsAsFactors = FALSE)
cscg.archaea=read.table(paste(cscg.dir,'cscg_archaea.tab', sep='/'),header=T,  stringsAsFactors = FALSE)
cscg=data.frame(pfam=c(cscg.bacteria$PFAM_ID,cscg.archaea$PFAM_ID),
                domain=c(rep('bacteria',nrow(cscg.bacteria)),rep('archaea', nrow(cscg.archaea))),
                duplicated=rep(0,nrow(cscg.bacteria)+(nrow(cscg.archaea))),
                absent=rep(0,nrow(cscg.bacteria)+(nrow(cscg.archaea))),
                duplication_rate=rep(0,nrow(cscg.bacteria)+nrow(cscg.archaea)),
                absence_rate=rep(0,nrow(cscg.bacteria)+nrow(cscg.archaea)),
                stringsAsFactors=FALSE)

progenomes.dir='/mnt/nfsdata/dev/florian/download-progenomes/progenomes'
progenomes.projects=basename(list.dirs(progenomes.dir))[-1]
progenomes.num_projects=length(progenomes.projects)

progenomes.cscg =data.frame(domain=rep('',progenomes.num_projects),
                 num_found=rep(0,progenomes.num_projects),
                 num_unique=rep(0,progenomes.num_projects),
                 num_duplicated=rep(0,progenomes.num_projects),
                 perc_found=rep(0,progenomes.num_projects),
                 stringsAsFactors=FALSE) 

rownames(progenomes.cscg)=progenomes.projects

for (project in rownames(progenomes.cscg))
{
  project.cscg.bacteria=read.table(paste(progenomes.dir, project, 'proteins.bacteria_cscg.tab',sep='/'),header=T, stringsAsFactors = FALSE)
  project.cscg.archaea=read.table(paste(progenomes.dir, project, 'proteins.archaea_cscg.tab',sep='/'),header=T, stringsAsFactors = FALSE)
  
  project.num_cscg.bacteria=length(unique(project.cscg.bacteria$pfam_name))
  project.num_cscg.archaea=length(unique(project.cscg.archaea$pfam_name))
  
  if (project.num_cscg.archaea > project.num_cscg.bacteria)
  {
    progenomes.cscg[project,]$domain="archaea"
    progenomes.cscg[project,]$num_found=project.num_cscg.archaea
    progenomes.cscg[project,]$num_duplicated=sum(duplicated(project.cscg.archaea$pfam_name))
    progenomes.cscg[project,]$num_unique=progenomes.cscg[project,]$num_found-progenomes.cscg[project,]$num_duplicated
    progenomes.cscg[project,]$perc_found=100*progenomes.cscg[project,]$num_found/sum(cscg$domain=='archaea')
    
    if (progenomes.cscg[project,]$perc_found >= 0.9)
    {
      for (duplicated_cscg in project.cscg.archaea$pfam_name[duplicated(project.cscg.archaea$pfam_name)])
      {
        cscg$duplicated[cscg$domain=='archaea' & cscg$pfam==duplicated_cscg]=
          cscg$duplicated[cscg$domain=='archaea' & cscg$pfam==duplicated_cscg]+1
      }
      
      for (absent_cscg in cscg$pfam[!(cscg$pfam %in% project.cscg.archaea$pfam_name) & cscg$domain == 'archaea'])
      {
        cscg$absent[cscg$domain=='archaea' & cscg$pfam==absent_cscg]=
          cscg$absent[cscg$domain=='archaea' & cscg$pfam==absent_cscg]+1
      }
    }
  } else {
    progenomes.cscg[project,]$domain="bacteria"
    progenomes.cscg[project,]$num_found=project.num_cscg.bacteria
    progenomes.cscg[project,]$num_duplicated=sum(duplicated(project.cscg.bacteria$pfam_name))
    progenomes.cscg[project,]$num_unique=progenomes.cscg[project,]$num_found-progenomes.cscg[project,]$num_duplicated
    progenomes.cscg[project,]$perc_found=100*progenomes.cscg[project,]$num_found/sum(cscg$domain=='bacteria')
    
    if (progenomes.cscg[project,]$perc_found  >= 0.9)
    {
      for (duplicated_cscg in project.cscg.bacteria$pfam_name[duplicated(project.cscg.bacteria$pfam_name)])
      {
         cscg$duplicated[cscg$domain=='bacteria' & cscg$pfam==duplicated_cscg]=
           cscg$duplicated[cscg$domain=='bacteria' & cscg$pfam==duplicated_cscg]+1
      }
      
      for (absent_cscg in cscg$pfam[!(cscg$pfam %in% project.cscg.bacteria$pfam_name) & cscg$domain == 'bacteria'])
      {
        cscg$absent[cscg$domain=='bacteria' & cscg$pfam==absent_cscg]=
          cscg$absent[cscg$domain=='bacteria' & cscg$pfam==absent_cscg]+1
      }
    }
  }
}

cscg$duplication_rate[cscg$domain=='archaea']=
  100*cscg$duplicated[cscg$domain=='archaea']/sum(progenomes.cscg$perc_found >= 0.9 & progenomes.cscg$domain=='archaea')

cscg$duplication_rate[cscg$domain=='bacteria']=
  100*cscg$duplicated[cscg$domain=='bacteria']/sum(progenomes.cscg$perc_found >= 0.9 & progenomes.cscg$domain=='bacteria')

cscg$absence_rate[cscg$domain=='archaea']=
  100*cscg$absent[cscg$domain=='archaea']/sum(progenomes.cscg$perc_found >= 0.9 & progenomes.cscg$domain=='archaea')

cscg$absence_rate[cscg$domain=='bacteria']=
  100*cscg$absent[cscg$domain=='bacteria']/sum(progenomes.cscg$perc_found >= 0.9 & progenomes.cscg$domain=='bacteria')

ggplot(progenomes.cscg, aes(factor(domain), fill = factor(domain), perc_found)) +
  geom_boxplot() +
  xlab("Domain") +
  ylab("% of Conserved Single Copy Genes found") +
  theme_bw() +
  theme(legend.position="none",
        text = element_text(size = 20)) 

ggplot(cscg, aes(x=100-duplication_rate, y=100-absence_rate, shape = domain)) +
  geom_point(aes(colour = domain), size = 3) +
  theme_bw() +
  theme(text = element_text(size = 20)) 