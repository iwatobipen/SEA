import math
import pytest
from sea.chemistry import normalize, raw_score, molecule, InputError
from sea.statistics import original_tail, significance, calibrate


def test_normalization_and_similarity():
    mols,rows=normalize(['CCO','OCC','CCO.[Na+]'])
    assert len(mols)==1 and sum(r['duplicate'] for r in rows)==2
    left,_=normalize(['CCO','CCCO'])
    right,_=normalize(['CCCCO','CCO'])
    a=raw_score(left,right,.57)
    b=raw_score(right,left,.57)
    assert a[0]==pytest.approx(b[0])
    assert a[1]>=1 and a[2][0]['tanimoto']==1
    assert raw_score([molecule('CCO')],[molecule('CCO')],1)[0]==1
    assert raw_score([molecule('CCO')],[molecule('c1ccccc1')],1)[0]==0
    with pytest.raises(InputError): normalize([])
    with pytest.raises(InputError) as exc: normalize(['CCO','wrong'])
    assert exc.value.args[0][0]['row']==2


def test_original_equations_precision():
    # Table 5 constants; Supplementary Methods uses rounded arithmetic in its example.
    z=(772.25-.000424*54648)/(.00449*54648**.665)
    assert 117<z<120
    p,logp=original_tail(119)
    assert p==pytest.approx(2.92e-67,rel=.02)
    p,logp=original_tail(1000)
    assert p==0 and math.isfinite(logp)
    assert original_tail(-100)[0]==1


def test_empirical_reproducibility():
    pool=[molecule(s) for s in ['CO','CCO','CCCO','CCCCO','CCCCC','CCCCCC','CCN','CCCN','CCCl','CCBr']]
    model={'threshold':.57,'std_b':0,'std_c':.5,'mean_a':0,'evd_usable':False,'seed':42,'empirical_trials':99,'sizes':[1,10]}
    a=significance(2,1,2,pool,model,3)
    assert a==significance(2,1,2,pool,model,3)
    assert a['p_value']>=.01 and a['e_value']==a['p_value']*3
    assert a['method']=='monte_carlo' and a['z_score'] is None
    assert significance(0,1,2,pool,model)['p_value']==1
    with pytest.raises(ValueError): significance(1,1,20,pool,model)


def test_gumbel_path_preserves_log_tail():
    model={'threshold':.57,'std_b':1.,'std_c':.5,'mean_a':0.,'evd_usable':True,'seed':42,'empirical_trials':99,'sizes':[1,10,100], 'gumbel_loc':0.,'gumbel_scale':1.}
    result=significance(10000,10,10,[],model,20)
    assert result['method']=='gumbel'
    assert result['p_value']==0 and result['log_p_value']==-1000
    assert result['log_e_value']==pytest.approx(-1000+math.log(20))


def test_zero_score_is_exact_without_sampling(monkeypatch):
    import sea.statistics as statistics
    def unexpected(*args, **kwargs):
        raise AssertionError('Zero scores must not run Monte Carlo sampling')
    monkeypatch.setattr(statistics.np.random, 'default_rng', unexpected)
    model={'threshold':.57,'std_b':1.,'std_c':.5,'mean_a':0.,'evd_usable':False,'seed':42,'empirical_trials':999,'sizes':[1,10]}
    result=significance(0,1,10,[],model,1517)
    assert result['method']=='exact_zero'
    assert result['p_value']==1 and result['log_p_value']==0
    assert result['e_value']==1517 and result['p_resolution'] is None
    # The nonnegative support also applies when an EVD model is available.
    model['evd_usable']=True
    assert significance(0,10,10,[],model)['method']=='exact_zero'
