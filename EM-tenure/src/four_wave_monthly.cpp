#include <Rcpp.h>
#include <array>
#include <cmath>
#include <unordered_map>
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

namespace {
const double neginf = -std::numeric_limits<double>::infinity();
double ladd(double a, double b) {
  if (a == neginf) return b;
  if (b == neginf) return a;
  double m = std::max(a, b);
  return m + std::log1p(std::exp(std::min(a,b)-m));
}
double mix(double a, double b, double q) {
  if (q == 0) return a;
  if (q == 1) return b;
  return ladd(std::log1p(-q)+a, std::log(q)+b);
}
int mod12(int x) { return ((x % 12) + 12) % 12; }
struct Spell { double ll, changes; };
struct SegmentHash {
  size_t operator()(const std::array<int,8>& key) const {
    size_t h=0;
    for(int x:key) h^=std::hash<int>{}(x)+0x9e3779b9+(h<<6)+(h>>2);
    return h;
  }
};
}

// Probability tables are constructed by the reference R routines. This
// kernel only enumerates the unchanged histories, job resets and report types.
// [[Rcpp::export]]
List four_wave_monthly_kernel(List data, List tables, NumericVector p,
                             bool posterior = true, int threads = 4) {
  IntegerMatrix y=data["y"], g=data["month"], d=data["category"], im=data["interview"];
  NumericVector w=data["weight"], raw=tables["raw"], norm=tables["normalizers"],
    dm=tables["timegap_marginal"], exit=tables["exit"], entry=tables["entry"];
  NumericMatrix base=tables["calendar"], within=tables["within"],
    dp=tables["timegap_pair"], dj3=tables["timegap_joint3"], dj4=tables["timegap_joint4"];
  bool clock_joint=tables["timegap_clock_joint"];
  std::array<NumericMatrix,4> clocks={{tables["timegap_clock1"],tables["timegap_clock2"],
    tables["timegap_clock3"],tables["timegap_clock4"]}};
  double alpha=p["alpha"], pi=p["pi"], q=p["job_change"], heap=p["heaping"],
    year=p["year"], cleanrev=p["clean_revision"], exact=p["exact"], local=p["local"];
  std::array<double,2> eps={{p["eps_reliable"],p["eps_unreliable"]}};
  double exit_missing=tables["exit_missing"], entry_missing=tables["entry_missing"];
  int n=y.nrow();
  NumericVector rowll(n), reliability(posterior?n:0), jobnum(posterior?n:0), jobden(posterior?n:0);
  NumericMatrix gamma(posterior?n:0,16);
  double total=0;
  checkUserInterrupt();
  #pragma omp parallel num_threads(threads)
  {
  std::unordered_map<std::array<int,8>,double,SegmentHash> segment_cache[4];
  std::unordered_map<std::array<int,8>,Spell,SegmentHash> spell_cache[4];
  #pragma omp for schedule(static)
  for (int i=0;i<n;++i) {
    std::array<double,4> gross;
    for (int t=0;t<4;++t) {
      gross[t]=0;
      if (y(i,t)) {
        int m=g(i,t), j=im(i,t)-1;
        double heaped=(m%12==j)?raw[m]-norm[j]:neginf;
        gross[t]=mix(base(m,j),heaped,heap);
      }
    }
    double kernels[32], changes[32], opportunities[16];
    for (int cl=0;cl<2;++cl) {
      double seg[2][4][4];
      for (int mode=0;mode<2;++mode) for (int a=0;a<4;++a) for (int b=a;b<4;++b) {
        std::array<int,8> key; key.fill(-2);
        for(int t=a;t<=b;++t) {
          key[t-a]=y(i,t)?g(i,t):-1;
          key[4+t-a]=y(i,t)?im(i,t):0;
        }
        auto& cache=segment_cache[2*cl+mode];
        auto found=cache.find(key);
        if(found!=cache.end()) {seg[mode][a][b]=found->second;continue;}
        int obs[4], k=0;
        for(int t=a;t<=b;++t) if(y(i,t)) obs[k++]=t;
        if (!k) { seg[mode][a][b]=0; continue; }
        double gross_mass[4][2], clean_mass[4], type_weight[5];
        int anchor_month[4];
        for(int nx=0;nx<=k;++nx)
          type_weight[nx]=nx*std::log(eps[cl])+(k-nx)*std::log1p(-eps[cl]);
        for(int j=0;j<k;++j) {
            int t=obs[j];
            int m0=g(i,t)-3*(t-a);
            anchor_month[j]=m0;
            clean_mass[j]=m0<0?neginf:mode==1?
              (m0<=2?within(m0,im(i,t)-1):neginf):
              base(m0,mod12(im(i,t)-1-3*(t-a)));
            gross_mass[j][0]=gross_mass[j][1]=gross[t];
            if(j>0) {
              int prev=obs[j-1], gap=3*(t-prev), expected=g(i,prev)+gap;
              double yr=neginf, loc=neginf;
              if(year>0 || cleanrev>0) {
                int required=mod12(im(i,t)-im(i,prev)+g(i,prev));
                double removed=(expected%12==required)?
                  std::min(std::exp(raw[expected]-norm[required]),1-1e-15):0;
                yr=(g(i,t)%12==required && g(i,t)!=expected)?
                  raw[g(i,t)]-norm[required]-std::log1p(-removed):neginf;
              }
              if(local>0) {
                int delta=g(i,t)-expected;
                double denominator=0;
                for(int v=-6;v<=6;++v) if(v!=0 && expected+v>=0)
                  denominator+=std::exp(-std::abs(v)/3.0);
                loc=(delta!=0 && std::abs(delta)<=6)?
                  -std::abs(delta)/3.0-std::log(denominator):neginf;
              }
              for(int previous_type=0;previous_type<2;++previous_type) {
                double mass=mix(gross[t],yr,previous_type?year:cleanrev);
                mass=mix(mass,loc,local);
                gross_mass[j][previous_type]=mix(mass,g(i,t)==expected?0:neginf,exact);
              }
            }
        }
        double ll=neginf;
        for(int mask=0;mask<(1<<k);++mask) {
          int nx=0, anchor=-1;
          for(int j=0;j<k;++j) {
            if(mask & (1<<j)) ++nx; else if(anchor<0) anchor=j;
          }
          double term=type_weight[nx];
          for(int j=0;j<k;++j) if(mask & (1<<j))
            term+=gross_mass[j][j>0?((mask>>(j-1))&1):0];
          if(anchor>=0) {
            int m0=anchor_month[anchor];
            bool valid=m0>=0;
            for(int j=0;j<k;++j) if(!(mask & (1<<j)))
              valid=valid && g(i,obs[j])==m0+3*(obs[j]-a);
            term+=valid?clean_mass[anchor]:neginf;
          }
          ll=ladd(ll,term);
        }
        seg[mode][a][b]=ll;
        cache.emplace(key,ll);
      }
      Spell spells[4][4];
      for(int a=0;a<4;++a) for(int b=a;b<4;++b) {
        std::array<int,8> key; key.fill(-2);
        for(int t=a;t<=b;++t) {
          key[t-a]=y(i,t)?g(i,t):-1;
          key[4+t-a]=y(i,t)?im(i,t):0;
        }
        auto& cache=spell_cache[2*cl+(a>0)];
        auto found=cache.find(key);
        if(found!=cache.end()) {spells[a][b]=found->second;continue;}
        int L=b-a+1;
        double ll=neginf, change_mass=neginf;
        for(int mask=0;mask<(1<<(L-1));++mask) {
          int count=0;
          for(int j=0;j<L-1;++j) count+=(mask>>j)&1;
          if(q==0 && count) continue;
          double term=count?count*std::log(q):0;
          term+=(L-1-count)*std::log1p(-q);
          int start=a;
          for(int t=a;t<=b;++t) if(t==b || (mask & (1<<(t-a)))) {
            term+=seg[(start==a)?(a>0):1][start][t]; start=t+1;
          }
          ll=ladd(ll,term);
          if(count) change_mass=ladd(change_mass,term+std::log(double(count)));
        }
        spells[a][b]={ll,change_mass==neginf?0:std::exp(change_mass-ll)};
        cache.emplace(key,spells[a][b]);
      }
      for(int h=0;h<16;++h) {
        double value=(h&1)?std::log(alpha):std::log1p(-alpha), jc=0, jo=0;
        int mismatch=0;
        for(int t=0;t<4;++t) mismatch+=(((h>>t)&1)!=y(i,t));
        if(pi==0 && mismatch) {kernels[cl*16+h]=neginf;changes[cl*16+h]=0;opportunities[h]=0;continue;}
        value+=(mismatch?mismatch*std::log(pi):0)+(4-mismatch)*std::log1p(-pi);
        for(int t=0;t<3;++t) {
          bool employed=(h>>t)&1;
          double risk=employed?(y(i,t)?exit[g(i,t)]:exit_missing):
            (!y(i,t)?entry[d(i,t)-1]:entry_missing);
          bool change=employed!=bool((h>>(t+1))&1);
          value+=change?std::log(risk):std::log1p(-risk);
        }
        bool covered[4]={false,false,false,false};
        for(int a=0;a<4;) {
          int state=(h>>a)&1, b=a;
          while(b+1<4 && ((h>>(b+1))&1)==state) ++b;
          if(state) {value+=spells[a][b].ll;jc+=spells[a][b].changes;jo+=b-a;}
          else {
            if(clock_joint) {
              int index=0,power=1;
              for(int t=a;t<=b;++t) {
                value+=gross[t];
                index+=(y(i,t)?0:d(i,t))*power;
                power*=8;covered[t]=true;
              }
              value+=clocks[b-a](index,cl);
            } else {
            bool all=true;int index=0,power=1;
            for(int t=a;t<=b;++t) {
              value+=gross[t];all=all && !y(i,t);
              if(!y(i,t)) index+=(d(i,t)-1)*power;
              power*=7;
            }
            if(b-a+1>=3 && all) {
              value+=(b-a+1==3)?dj3(index,cl):dj4(index,cl);
              for(int t=a;t<=b;++t) covered[t]=true;
            }
            }
          }
          a=b+1;
        }
        for(int t=0;t<4;++t) if(!y(i,t) && !covered[t]) {
          bool conditional=t>0 && !((h>>(t-1))&1) && !((h>>t)&1) &&
            !y(i,t-1) && !covered[t-1];
          value+=conditional?dp((d(i,t-1)-1)+7*(d(i,t)-1),cl):dm[d(i,t)-1];
        }
        kernels[cl*16+h]=value-std::log(2.0);changes[cl*16+h]=jc;opportunities[h]=jo;
      }
    }
    double ll=neginf;
    for(int j=0;j<32;++j) ll=ladd(ll,kernels[j]);
    rowll[i]=ll;
    if(posterior && std::isfinite(ll)) for(int j=0;j<32;++j) {
      double prob=std::exp(kernels[j]-ll);
      gamma(i,j%16)+=prob; reliability[i]+=(j>=16)?prob:0;
      jobnum[i]+=prob*changes[j];jobden[i]+=prob*opportunities[j%16];
    }
  }
  }
  checkUserInterrupt();
  for(int i=0;i<n;++i) total+=w[i]*rowll[i];
  return List::create(_["loglik"]=total,_["row_loglik"]=rowll,_["gamma"]=gamma,
    _["duration_reliability_posterior"]=reliability,
    _["job_change_posterior"]=List::create(_["expected_changes"]=jobnum,_["opportunities"]=jobden));
}
